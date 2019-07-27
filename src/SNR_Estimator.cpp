#include "SNR_Estimator.h"

SNR_Estimator::SNR_Estimator() {
    norm2_ = 0.0;
    norm4_ = 0.0;
    sampleCount_ = 0;
    lastEstimate_ = 0.0;
}

SNR_Estimator::~SNR_Estimator() {

}

void SNR_Estimator::input(std::complex<double> carrier) {
    double norm2 = carrier.real() * carrier.real() + carrier.imag() * carrier.imag();
    norm2_ += norm2;
    norm4_ += norm2 * norm2;
    ++sampleCount_;
}

double SNR_Estimator::estimate(void) {
    if (sampleCount_ < 500)
        return lastEstimate_;
    norm2_ /= sampleCount_;
    norm4_ /= sampleCount_;
    double powerD = sqrt(2.0 * norm2_ * norm2_ - norm4_);
    double powerN = norm2_ - powerD;
    norm2_ = 0.0;
    norm4_ = 0.0;
    sampleCount_ = 0;
    lastEstimate_ = powerD / powerN;
    if (!std::isfinite(lastEstimate_))
        lastEstimate_ = 0.0;
    return lastEstimate_;
}
