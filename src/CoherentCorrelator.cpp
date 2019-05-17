#include "CoherentCorrelator.h"

CoherentCorrelator::CoherentCorrelator(double fs, unsigned prn, bool recordOn) : caCode_(prn), fs_(fs) {
	recorderOn_ = recordOn;
	if (recorderOn_)
		recorder_.open("correlator.bin", std::ofstream::binary);
	ready_ = false;
	integrationCount_ = 0;
	lastCaChipNum_ = 0;
}

CoherentCorrelator::~CoherentCorrelator() {
	if (recorderOn_)
		recorder_.close();
}

double CoherentCorrelator::integrate(fftwVector &data, double codeFreqOffset, double timeOffset, unsigned integrationCount, bool *ready){
	timeOffset += 1e-3;

	for (unsigned i = 0; i < data.size(); ++i){
		double time = i / fs_ + timeOffset;
		int chipNumber = static_cast<int>(round(time * (CHIP_RATE + codeFreqOffset))) % CA_CODE_LENGTH;

		if (lastCaChipNum_ == 1022 && chipNumber == 0){
			timeOffset -= 1e-3;
			++integrationCount_;
		}
		lastCaChipNum_ = chipNumber;

		if (integrationCount_ == integrationCount){
			ready_ = true;
			readyVal_ = currentIntegration_;
			currentIntegration_ = std::complex<double>(0.0,0.0);
			integrationCount_ = 0;
		}

		int chip = caCode_[chipNumber] == 1 ? 1 : -1;
		currentIntegration_ += data[i] * std::complex<double>(chip, 0.0);
		std::complex<double> rec(chip, 0.0);
	}

	if (ready)
		*ready &= ready_;
	return timeOffset;
}

std::complex<double> CoherentCorrelator::dump(void){
	if (ready_){
		ready_ = false;
		return readyVal_;
	}
	std::complex<double> retval = currentIntegration_;
	currentIntegration_ = std::complex<double>(0.0,0.0);
	integrationCount_ = 0;
	return retval;
}

std::complex<double> CoherentCorrelator::val(void){
	return currentIntegration_;
}

double CoherentCorrelator::integrationTime(void){
	return integrationCount_ * 1e-3;
}
