#include "SignalTracker.h"
#include <unistd.h>

const double CARRIER_CODE_SF = 1.023/1575.42;

SignalTracker::SignalTracker(double fs, unsigned prn, SearchResult searchResult, double searchFreqOffset) : runThread_(true),
early_(fs, prn),
prompt_(fs, prn, true),
late_(fs, prn),
carrierCorrelator_(fs, prn),
carrierPhaseLPF_(0.0247),
carrierFreqLPF_(0.0247),
processingThread_(&SignalTracker::threadFunction, this) {
	prn_ = prn;
	fs_ = fs;

	timeSinceStart_ = 0.0;
	processedPackets_ = 0;

	carrierPhase_        = 0.0;
	carrierFreq_         = searchResult.baseBandFreq + searchFreqOffset;
	codetime_            = searchResult.sampleOffset / fs_;
	codeFreq_            = (searchResult.baseBandFreq + searchFreqOffset) * CARRIER_CODE_SF;
	integrationLength_   = 1;

	lastCarrier_ = complex<double>(1.0, 0.0);

	state_ = closingCarrierFLL;

	codeRecorder_.open("code" + std::to_string(prn) + '_' + std::to_string((int)searchFreqOffset) + ".bin", std::ofstream::binary);
	carrierRecorder_.open("carrier" + std::to_string(prn) + '_' + std::to_string((int)searchFreqOffset) + ".bin", std::ofstream::binary);
}

SignalTracker::~SignalTracker() {
	runThread_ = false;
	processingThread_.join();
	codeRecorder_.close();
	carrierRecorder_.close();
}

State SignalTracker::state(void){
	return state_;
}

double SignalTracker::CNoEst(void){
	return snrEstimator_.estimate() / (CA_CODE_TIME * integrationLength_);
}

void SignalTracker::threadFunction(void){
	while (runThread_){
		fftwVector trackingData;
		trackingDataAccess_.lock();
		size_t size = trackingData_.size();
		trackingDataAccess_.unlock();
		if (size){
			trackingDataAccess_.lock();
			trackingData = std::move(trackingData_.front());
			trackingData_.pop_front();
			trackingDataAccess_.unlock();
		}
		else {
			usleep(1e3);
			continue;
		}

		frequencyShift(trackingData, -carrierFreq_);

		carrierCorrelator_.integrate(trackingData, codeFreq_, codetime_, 1);
		std::complex<double> carrier = carrierCorrelator_.dump();
		snrEstimator_.input(carrier);
		double CN0 = snrEstimator_.estimate() / (CA_CODE_TIME);

		complex<double> carrierError = carrier / lastCarrier_;
		double carrierErrorPhase = arg(carrierError);
		double carrierFreqHz = carrierErrorPhase / (CA_CODE_TIME * 2.0 * M_PI);

		complex<double> costasLoopError = carrier.real() > 0.0 ? carrier : carrier * complex<double>(-1.0,0.0);
		double costasLoopErrorPhase = arg(costasLoopError);
		double costasLoopFreqError  = costasLoopErrorPhase / (CA_CODE_TIME * 2.0 * M_PI);

		double costasLoopPhaseCorrection = state_ == closingCarrierFLL ? 0.0 : -costasLoopErrorPhase / 5.0;
		double fllFreqError = state_ == closingCarrierFLL ? carrierFreqHz : costasLoopFreqError;
		double fllGain = state_ == closingCarrierFLL ? .01 : .001;

		lastCarrier_ = carrier;
		carrierFllLock_.error(carrierError);
		carrierPhaseLPF_.iterate(costasLoopErrorPhase);
		carrierFreqLPF_.iterate(fllFreqError);

		switch (state_){
		case lossOfLock:
			break;
		case closingCarrierFLL:
			if (!codeLock_.optimisticLock())
				state_ = lossOfLock;
			if (carrierFllLock_.pessimisticLock())
				state_ = closingCarrierPLL;
			break;
		case closingCarrierPLL:
			integrationLength_ = 1;
			carrierPllLock_.error(costasLoopError);
			if (!carrierFllLock_.optimisticLock())
				state_ = closingCarrierFLL;
			if (carrierPllLock_.pessimisticLock())
				state_ = findingNavBitEdge;
			break;
		case findingNavBitEdge:
			carrierPllLock_.error(costasLoopError);
			if (!carrierPllLock_.optimisticLock())
				state_ = closingCarrierPLL;
			if (fabs(carrierErrorPhase) > 0.5 * M_PI){
				if (navBitEdgeDetector_.push_back(processedPackets_)){
					integrationLength_ = 20;
					state_ = fullTrack;
				}
			}
			break;
		case fullTrack:
			break;
		}

		carrierFreq_  += fllGain * fllFreqError;
		carrierPhase_ += costasLoopPhaseCorrection;
		//Aid code tracking
		codeFreq_     -= fllGain * fllFreqError * CARRIER_CODE_SF;
		codetime_     -= costasLoopPhaseCorrection / (2.0 * M_PI * 1.57542e9) * CARRIER_CODE_SF;
		lastCarrier_  *= complex<double>(cos(costasLoopPhaseCorrection), sin(costasLoopPhaseCorrection));

		double state = state_;
		double meanCarrierPhaseErr = carrierPhaseLPF_.last();
		double meanCarrierFreqErr = carrierFreqLPF_.last();
		carrierRecorder_.write((char*)&timeSinceStart_, sizeof(timeSinceStart_));
		carrierRecorder_.write((char*)&carrierError, sizeof(carrier));
		carrierRecorder_.write((char*)&costasLoopErrorPhase, sizeof(costasLoopErrorPhase));
		carrierRecorder_.write((char*)&meanCarrierPhaseErr, sizeof(meanCarrierPhaseErr));
		carrierRecorder_.write((char*)&fllFreqError, sizeof(fllFreqError));
		carrierRecorder_.write((char*)&meanCarrierFreqErr, sizeof(meanCarrierFreqErr));
		carrierRecorder_.write((char*)&carrierFreq_, sizeof(carrierFreq_));
		carrierRecorder_.write((char*)&CN0, sizeof(CN0));
		carrierRecorder_.write((char*)&state, sizeof(state));
		double ilast = carrierFllLock_.lastInPhase();
		double qlast = carrierFllLock_.lastQuadPhase();
		double p1    = carrierFllLock_.pCount1();
		double p2    = carrierFllLock_.pCount2();
		double opt   = carrierFllLock_.optimisticLock();
		double pes   = carrierFllLock_.pessimisticLock();
		carrierRecorder_.write((char*)&ilast, sizeof(ilast));
		carrierRecorder_.write((char*)&qlast, sizeof(qlast));
		carrierRecorder_.write((char*)&p1, sizeof(p1));
		carrierRecorder_.write((char*)&p2, sizeof(p2));
		carrierRecorder_.write((char*)&opt, sizeof(opt));
		carrierRecorder_.write((char*)&pes, sizeof(pes));

		bool ready = true;
		early_.integrate(trackingData, codeFreq_, codetime_ - CORR_OFFSET, integrationLength_);
		codetime_ = prompt_.integrate(trackingData, codeFreq_, codetime_, integrationLength_, &ready);
		late_.integrate(trackingData, codeFreq_, codetime_ + CORR_OFFSET, integrationLength_);

		if (ready) {
			std::complex<double> early = early_.dump();
			std::complex<double> prompt = prompt_.dump();
			std::complex<double> late = late_.dump();
			double magEarly  = abs(early);
			double magPrompt = abs(prompt);
			double magLate   = abs(late);
			double codeErrorChips = (magEarly - magLate) / (2.0 * (magEarly + magLate));
			double codeErrorTime = codeErrorChips / CHIP_RATE;
			codeLock_.error(magEarly, magPrompt, magLate);

			codetime_ -= 0.04 * codeErrorTime;
			codeFreq_ -= 0.04 * codeErrorChips;

			double lpEarly = codeLock_.early();
			double lpPrompt = codeLock_.prompt();
			double lpLate = codeLock_.late();
			codeRecorder_.write((char*)&timeSinceStart_, sizeof(timeSinceStart_));
			codeRecorder_.write((char*)&early, sizeof(early));
			codeRecorder_.write((char*)&prompt, sizeof(prompt));
			codeRecorder_.write((char*)&late, sizeof(late));
			codeRecorder_.write((char*)&codetime_, sizeof(codetime_));
			codeRecorder_.write((char*)&codeFreq_, sizeof(codeFreq_));
			codeRecorder_.write((char*)&lpEarly, sizeof(lpEarly));
			codeRecorder_.write((char*)&lpPrompt, sizeof(lpPrompt));
			codeRecorder_.write((char*)&lpLate, sizeof(lpLate));

			if (integrationLength_ == 20){
				bool flipBits;
				lnav_data_.navBit(2 * (fabs(arg(prompt)) < (M_PI / 2.0)) - 1, flipBits);
				if (flipBits)
					carrierPhase_ += M_PI;
			}
		}

		codetime_ += codeFreq_ * CA_CODE_TIME * CHIP_TIME;
		timeSinceStart_ += trackingData.size() / fs_;
		++processedPackets_;
	}
}

unsigned SignalTracker::prn(void){
	return prn_;
}

double SignalTracker::transmitTime(void){
	double navTime = lnav_data_.timeOfLastNavBit();
	if (navTime == -1.0)
		return -1.0;
	return navTime + prompt_.integrationTime() + codetime_; //CA_CODE_TIME -
}

Vector3d SignalTracker::satellitePosition(double timeOfWeek){
	if (timeOfWeek == -1.0)
		return Vector3d(0.0, 0.0, 0.0);
	return lnav_data_.satellitePosition(timeOfWeek);
}

complex<double> SignalTracker::latLong(double timeOfWeek){
	Vector3d ecef = lnav_data_.satellitePosition(timeOfWeek);
	complex<double> retVal;
	retVal.real(atan2(ecef.z(), sqrt(ecef.x() * ecef.x() + ecef.y() * ecef.y())) * 180.0 / M_PI);
	retVal.imag(atan2(ecef.y(), ecef.x()) * 180.0 / M_PI);
	return retVal;
}

bool SignalTracker::processSamples(fftwVector trackingData){
	trackingDataAccess_.lock();
	trackingData_.push_back(std::move(trackingData));
	trackingDataAccess_.unlock();
	return state_ != lossOfLock;
}

void SignalTracker::sync(void){
	while(1){
		trackingDataAccess_.lock();
		if (trackingData_.size()){
			trackingDataAccess_.unlock();
			usleep(1e3);
		}
		else {
			trackingDataAccess_.unlock();
			break;
		}
	}
}

void SignalTracker::frequencyShift(fftwVector &data, double frequencyShift_Hz){
	carrierPhase_ = fmod(carrierPhase_, 2.0 * M_PI);
	for (unsigned i = 0; i < data.size(); ++i){
		carrierPhase_ += 2.0 * M_PI * frequencyShift_Hz / fs_;
		std::complex<double> mult(cos(carrierPhase_), sin(carrierPhase_));
		data[i] *= mult;
	}
}
