#include "CodeLockDetector.h"

const double LPF_FCUTOFF = .0039805;
const unsigned OPTIMISTIC_THRESH = 240;

CodeLockDetector::CodeLockDetector() : earlyLPF_(LPF_FCUTOFF), promptLPF_(LPF_FCUTOFF), lateLPF_(LPF_FCUTOFF) {
    optimisticCount_ = 0;
}

CodeLockDetector::~CodeLockDetector() {

}

void CodeLockDetector::error(double magEarly, double magPrompt, double magLate) {
    double early = earlyLPF_.iterate(magEarly);
    double prompt = promptLPF_.iterate(magPrompt);
    double late = lateLPF_.iterate(magLate);

    if (early > prompt || late > prompt)
        ++optimisticCount_;
    else
        optimisticCount_ = 0;
}

bool CodeLockDetector::optimisticLock(void) {
    return optimisticCount_ < OPTIMISTIC_THRESH;
}
