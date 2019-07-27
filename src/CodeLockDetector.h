#ifndef SRC_CODELOCKDETECTOR_H_
#define SRC_CODELOCKDETECTOR_H_

#include "LowPassFilter.h"

/*
 * This is a paired-down version of the lock detector, specifically for detecting code-lock. Since code-lock is assumed at the start of
 * the signal tracker, and a loss of code-lock is considered non-recoverable for the signal tracker, this is only a loss-of-lock detector,
 * and hence no pessimistic side for detecting start of lock.
 */

class CodeLockDetector {
public:
    CodeLockDetector();

    virtual ~CodeLockDetector();

    void error(double magEarly, double magPrompt, double magLate);

    bool optimisticLock(void);

    double early(void) const {
        return earlyLPF_.last();
    }

    double prompt(void) const {
        return promptLPF_.last();
    }

    double late(void) const {
        return lateLPF_.last();
    }

private:
    LowPassFilter earlyLPF_;
    LowPassFilter promptLPF_;
    LowPassFilter lateLPF_;

    unsigned optimisticCount_;
};

#endif /* SRC_CODELOCKDETECTOR_H_ */
