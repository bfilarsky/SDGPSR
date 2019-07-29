#include "LowPassFilter.h"
#include <cmath>
#include <stdexcept>

LowPassFilter::LowPassFilter(double fCutoff) {
    if (fCutoff <= 0.0 || fCutoff > 0.5)
        throw std::range_error("Cutoff frequency must be in the range (0,0.5]");
    k_ = 1.0 - exp(-2.0 * M_PI * fCutoff);
    last_ = 0.0;
}

LowPassFilter::~LowPassFilter() {

}

double LowPassFilter::iterate(double input) {
    last_ += k_ * (input - last_);
    return last_;
}

double LowPassFilter::last(void) const {
    return last_;
}
