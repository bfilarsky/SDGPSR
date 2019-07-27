#include "LockDetector.h"

const double LPF_FCUTOFF = .0039805;
const double K2 = 1.5;    //In phase divisor
const unsigned PESSIMISTIC_THRESH = 50;
const unsigned OPTIMISTIC_THRESH = 240;

LockDetector::LockDetector() : inPhaseLPF_(LPF_FCUTOFF), quadPhaseLPF_(LPF_FCUTOFF) {
    pessimisticCount_ = 0;
    optimisticCount_ = 0;
    optimisticLock_ = false;
    pessimisticLock_ = false;
}

LockDetector::~LockDetector() {

}

void LockDetector::error(std::complex<double> error) {
    double A = inPhaseLPF_.iterate(fabs(error.real())) / K2;
    double B = quadPhaseLPF_.iterate(fabs(error.imag()));
    if (A > B) {
        optimisticLock_ = true;
        optimisticCount_ = 0;
        if (pessimisticCount_ > PESSIMISTIC_THRESH)
            pessimisticLock_ = true;
        else
            ++pessimisticCount_;
    } else {
        pessimisticLock_ = false;
        pessimisticCount_ = 0;
        if (optimisticCount_ > OPTIMISTIC_THRESH)
            optimisticLock_ = false;
        else
            ++optimisticCount_;
    }
}

bool LockDetector::optimisticLock(void) const {
    return optimisticLock_;
}

bool LockDetector::pessimisticLock(void) const {
    return pessimisticLock_;
}

double LockDetector::lastInPhase(void) const {
    return inPhaseLPF_.last();
}

double LockDetector::lastQuadPhase(void) const {
    return quadPhaseLPF_.last();
}

unsigned LockDetector::pessimisticCount(void) const {
    return pessimisticCount_;
}

unsigned LockDetector::optimisticCount(void) const {
    return optimisticCount_;
}
