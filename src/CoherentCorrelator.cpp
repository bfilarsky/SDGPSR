#include "CoherentCorrelator.h"

CoherentCorrelator::CoherentCorrelator(double fs, unsigned prn) : caCode_(prn), fs_(fs) {
    ready_ = false;
    integrationCount_ = 0;
    lastCaChipNum_ = 0;
}

CoherentCorrelator::~CoherentCorrelator() {

}

double CoherentCorrelator::integrate(fftwVector &data, double codeFreqOffset, double timeOffset, unsigned integrationCount, bool *ready) {
    timeOffset += 1e-3;

    for (unsigned i = 0; i < data.size(); ++i) {
        double time = i / fs_ + timeOffset;
        int chipNumber = static_cast<int>(round(time * (CHIP_RATE + codeFreqOffset))) % CA_CODE_LENGTH;

        if (lastCaChipNum_ == 1022 && chipNumber == 0) {
            timeOffset -= 1e-3;
            ++integrationCount_;
        }
        lastCaChipNum_ = chipNumber;

        if (integrationCount_ == integrationCount) {
            ready_ = true;
            readyVal_ = currentIntegration_;
            currentIntegration_ = std::complex<double>(0.0, 0.0);
            integrationCount_ = 0;
        }

        currentIntegration_ += data[i] * std::complex<double>(caCode_[chipNumber], 0.0);
    }

    if (ready)
        *ready &= ready_;
    return timeOffset;
}

std::complex<double> CoherentCorrelator::dump(void) {
    if (ready_) {
        ready_ = false;
        return readyVal_;
    }
    std::complex<double> retval = currentIntegration_;
    currentIntegration_ = std::complex<double>(0.0, 0.0);
    integrationCount_ = 0;
    return retval;
}

std::complex<double> CoherentCorrelator::val(void) const {
    return currentIntegration_;
}

double CoherentCorrelator::integrationTime(void) const {
    return integrationCount_ * 1e-3;
}
