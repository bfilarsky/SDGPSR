#include "SignalTracker.h"
#include <chrono>
#include <thread>

const double CARRIER_CODE_SF  = 1.023 / 1575.42;
const double LPF_FCUTOFF      = .0039805;
const double FLL_GAIN_PULLIN  = .01;
const double FLL_GAIN_LOCKED  = .001;
const double CODE_PLL_P_GAIN  = .04;
const double CODE_PLL_I_GAIN  = .002;
const double COSTAS_LOOP_GAIN = 0.2;

SignalTracker::SignalTracker(double fs, unsigned prn, SearchResult searchResult) :
                                             runThread_(true),
                                             synced_(false),
                                             early_(fs, prn),
                                             prompt_(fs, prn),
                                             late_(fs, prn),
                                             carrierCorrelator_(fs, prn),
                                             navBitEdgeDetector_(5),
                                             carrierPhaseLPF_(LPF_FCUTOFF),
                                             carrierFreqLPF_(LPF_FCUTOFF),
                                             processingThread_(&SignalTracker::threadFunction, this) {
    prn_ = prn;
    fs_ = fs;

    timeSinceStart_ = 0.0;
    processedPackets_ = 0;

    carrierPhase_ = 0.0;
    carrierFreq_ = searchResult.baseBandFreq;
    codetime_ = searchResult.sampleOffset / fs_;
    codeFreq_ = (searchResult.baseBandFreq) * CARRIER_CODE_SF;
    integrationLength_ = 1;

    lastCarrier_ = std::complex<double>(1.0, 0.0);

    state_ = closingCarrierFLL;

#ifdef DEBUG_FILES
    codeRecorder_.open("code" + std::to_string(prn) +".bin", std::ofstream::binary);
    carrierRecorder_.open("carrier" + std::to_string(prn) + ".bin", std::ofstream::binary);
#endif
}

SignalTracker::~SignalTracker() {
    runThread_ = false;
    processingThread_.join();
}

State SignalTracker::state(void) {
    return state_;
}

double SignalTracker::CNoEst(void) {
    std::lock_guard<std::mutex> lock(snrMutex_);
    return snrEstimator_.estimate() / (CA_CODE_TIME * integrationLength_);
}

void SignalTracker::threadFunction(void) {
    while (runThread_) {
        fftwVector trackingData;
        trackingDataAccess_.lock();
        size_t size = trackingData_.size();
        trackingDataAccess_.unlock();
        if (size) {
            synced_ = false;
            trackingDataAccess_.lock();
            trackingData = std::move(trackingData_.front());
            trackingData_.pop_front();
            trackingDataAccess_.unlock();
        } else {
            synced_ = true;
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            continue;
        }

        frequencyShift(trackingData, -carrierFreq_);

        carrierCorrelator_.integrate(trackingData, codeFreq_, codetime_, 1);
        std::complex<double> carrier = carrierCorrelator_.dump();
        snrMutex_.lock();
        snrEstimator_.input(carrier);
        snrMutex_.unlock();

        switch (state_) {
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
                if (!carrierFllLock_.optimisticLock())
                    state_ = closingCarrierFLL;
                if (carrierPllLock_.pessimisticLock())
                    state_ = findingNavBitEdge;
                break;
            case findingNavBitEdge:
                if (!carrierPllLock_.optimisticLock())
                    state_ = closingCarrierPLL;
                //Since the PLL is locked, most of the energy should be in the in-phase
                //If a 180* flip has happened from a nav bit, carrier.real and lastCarrier_.real
                //will have opposite signs
                if (carrier.real() * lastCarrier_.real() < 0.0) {
                    if (navBitEdgeDetector_.push_back(processedPackets_)) {
                        integrationLength_ = 20;
                        state_ = fullTrack;
                    }
                }
                break;
                //TODO: Handle drops from below states
            case fullTrack:
                if (lnav_data_.valid())
                    state_ = fullNav;
                break;
            case fullNav:
                if (!lnav_data_.valid())
                    state_ = fullTrack;
                break;
        }

        std::complex<double> carrierRotation = carrier / lastCarrier_;
        double carrierErrorFreqHz = arg(carrierRotation) / (CA_CODE_TIME * 2.0 * M_PI);

        std::complex<double> costasLoopError = carrier.real() > 0.0 ? carrier : carrier * std::complex<double>(-1.0, 0.0);
        double costasLoopErrorPhase = arg(costasLoopError);
        double costasLoopFreqError = costasLoopErrorPhase / (CA_CODE_TIME * 2.0 * M_PI);

        lastCarrier_ = carrier;
        carrierFllLock_.error(carrierRotation);
        carrierPhaseLPF_.iterate(costasLoopErrorPhase);

        if (state_ == closingCarrierFLL){
            //Aid code tracking
            carrierFreqLPF_.iterate(carrierErrorFreqHz);
            carrierFreq_ += FLL_GAIN_PULLIN * carrierErrorFreqHz;
            codeFreq_    += FLL_GAIN_PULLIN * carrierErrorFreqHz * CARRIER_CODE_SF;
        }
        else {
            carrierPllLock_.error(costasLoopError);
            double costasLoopPhaseCorrection = -costasLoopErrorPhase * COSTAS_LOOP_GAIN;
            carrierPhase_ += costasLoopPhaseCorrection;
            codetime_     += costasLoopPhaseCorrection / (2.0 * M_PI * 1.57542e9) * CARRIER_CODE_SF;
            lastCarrier_  *= std::complex<double>(cos(costasLoopPhaseCorrection), sin(costasLoopPhaseCorrection));

            //Aid code tracking
            carrierFreqLPF_.iterate(costasLoopFreqError);
            carrierFreq_ += FLL_GAIN_LOCKED * costasLoopFreqError;
            codeFreq_    += FLL_GAIN_LOCKED * costasLoopFreqError * CARRIER_CODE_SF;
        }
        codetime_ += codeFreq_ * CA_CODE_TIME * CHIP_TIME;

#ifdef DEBUG_FILES
        double state = state_;
        double meanCarrierPhaseErr = carrierPhaseLPF_.last();
        double meanCarrierFreqErr = carrierFreqLPF_.last();
        double CN0 = snrEstimator_.estimate() / (CA_CODE_TIME);
        carrierRecorder_.write((char*) &timeSinceStart_, sizeof(timeSinceStart_));
        carrierRecorder_.write((char*) &carrierRotation, sizeof(carrier));
        carrierRecorder_.write((char*) &costasLoopErrorPhase, sizeof(costasLoopErrorPhase));
        carrierRecorder_.write((char*) &meanCarrierPhaseErr, sizeof(meanCarrierPhaseErr));
        carrierRecorder_.write((char*) &carrierErrorFreqHz, sizeof(carrierErrorFreqHz));
        carrierRecorder_.write((char*) &meanCarrierFreqErr, sizeof(meanCarrierFreqErr));
        carrierRecorder_.write((char*) &carrierFreq_, sizeof(carrierFreq_));
        carrierRecorder_.write((char*) &CN0, sizeof(CN0));
        carrierRecorder_.write((char*) &state, sizeof(state));
        double ilast = carrierFllLock_.lastInPhase();
        double qlast = carrierFllLock_.lastQuadPhase();
        double pesCount = carrierFllLock_.pessimisticCount();
        double optCount = carrierFllLock_.optimisticCount();
        double opt = carrierFllLock_.optimisticLock();
        double pes = carrierFllLock_.pessimisticLock();
        carrierRecorder_.write((char*) &ilast, sizeof(ilast));
        carrierRecorder_.write((char*) &qlast, sizeof(qlast));
        carrierRecorder_.write((char*) &pesCount, sizeof(pesCount));
        carrierRecorder_.write((char*) &optCount, sizeof(optCount));
        carrierRecorder_.write((char*) &opt, sizeof(opt));
        carrierRecorder_.write((char*) &pes, sizeof(pes));
#endif

        bool ready = true;
        early_.integrate(trackingData, codeFreq_, codetime_ - CORR_OFFSET, integrationLength_);
        codetime_ = prompt_.integrate(trackingData, codeFreq_, codetime_, integrationLength_, &ready);
        late_.integrate(trackingData, codeFreq_, codetime_ + CORR_OFFSET, integrationLength_);

        if (ready) {
            std::complex<double> early = early_.dump();
            std::complex<double> prompt = prompt_.dump();
            std::complex<double> late = late_.dump();
            double magEarly = abs(early);
            double magPrompt = abs(prompt);
            double magLate = abs(late);
            double codeErrorChips = (magLate - magEarly) / (2.0 * (magEarly + magLate));
            double codeErrorTime = codeErrorChips / CHIP_RATE;
            codeLock_.error(magEarly, magPrompt, magLate);

            codetime_ += CODE_PLL_P_GAIN * codeErrorTime;
            codeFreq_ += CODE_PLL_I_GAIN * codeErrorChips * integrationLength_;

#ifdef DEBUG_FILES
            double lpEarly = codeLock_.early();
            double lpPrompt = codeLock_.prompt();
            double lpLate = codeLock_.late();
            codeRecorder_.write((char*) &timeSinceStart_, sizeof(timeSinceStart_));
            codeRecorder_.write((char*) &early, sizeof(early));
            codeRecorder_.write((char*) &prompt, sizeof(prompt));
            codeRecorder_.write((char*) &late, sizeof(late));
            codeRecorder_.write((char*) &codetime_, sizeof(codetime_));
            codeRecorder_.write((char*) &codeFreq_, sizeof(codeFreq_));
            codeRecorder_.write((char*) &lpEarly, sizeof(lpEarly));
            codeRecorder_.write((char*) &lpPrompt, sizeof(lpPrompt));
            codeRecorder_.write((char*) &lpLate, sizeof(lpLate));
#endif

            if (state_ >= fullTrack) {
                bool flipBits;
                lnavMutex_.lock();
                lnav_data_.navBit(2 * (fabs(arg(prompt)) < (M_PI / 2.0)) - 1, flipBits);
                lnavMutex_.unlock();
                if (flipBits)
                    carrierPhase_ += M_PI;
            }
        }

        timeSinceStart_ += trackingData.size() / fs_;
        ++processedPackets_;
    }
}

unsigned SignalTracker::prn(void) {
    return prn_;
}

double SignalTracker::transmitTime(void) {
    lnavMutex_.lock();
    double navTime = lnav_data_.timeOfLastNavBit();
    lnavMutex_.unlock();
    if (navTime == -1.0)
        return -1.0;
    return navTime + prompt_.integrationTime() + codetime_;
}

Vector3d SignalTracker::satellitePosition(double timeOfWeek) {
    if (timeOfWeek == -1.0)
        return Vector3d(0.0, 0.0, 0.0);
    std::lock_guard<std::mutex> lock(lnavMutex_);
    return lnav_data_.satellitePosition(timeOfWeek);
}

State SignalTracker::processSamples(fftwVector trackingData) {
    trackingDataAccess_.lock();
    trackingData_.push_back(std::move(trackingData));
    trackingDataAccess_.unlock();
    return state_;
}

void SignalTracker::sync(void) {
    while (1) {
        trackingDataAccess_.lock();
        size_t dataSize = trackingData_.size();
        trackingDataAccess_.unlock();
        if (!synced_ || dataSize)
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        else
            return;
    }
}

void SignalTracker::frequencyShift(fftwVector &data, double frequencyShift_Hz) {
    carrierPhase_ = fmod(carrierPhase_, 2.0 * M_PI);
    for (unsigned i = 0; i < data.size(); ++i) {
        carrierPhase_ += 2.0 * M_PI * frequencyShift_Hz / fs_;
        std::complex<double> mult(cos(carrierPhase_), sin(carrierPhase_));
        data[i] *= mult;
    }
}
