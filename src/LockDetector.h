#ifndef SRC_LOCKDETECTOR_H_
#define SRC_LOCKDETECTOR_H_

#include <complex>
#include "LowPassFilter.h"

/*
 * Lock detector with optimistic and pessimistic decisions. This is taken from Kaplan: Understanding GPS - Principles and Applications 2nd edition, 5.11.
 * When the phase is locked, most of the error output should be on the in-phase component (no quad-phase component indicates no error).
 * This compares the low-pass filtered in-phase and quad-phase component - when the in-phase term / K2 is greater than the quad phase component, it
 * trends towards a locked state, otherwise it trends towards an unlocked state. The optimistic lock detection determines locked as soon as it trends toward lock,
 * while the pessimistic lock detection determines unlocked as soon as it trends towards an unlocked state. Conversely, the optimistic lock waits until OPTIMISTIC_THRESH
 * iterations have passed trending towards the unlocked direction before it gives up, while the pessimistic lock waits until PESSIMISTIC_THRESH before it determines
 * it is locked. Typically, the pessimistic lock is used to make a decision that the signal is locked, and the optimistic lock is used to make a decision on loss of lock.
 */

class LockDetector {
public:
    LockDetector();

    virtual ~LockDetector();

    //Add current error, where error is phase angle from 0
    void error(std::complex<double> error);

    //Return status of the optimistic lock - true for locked, false for unlocked
    bool optimisticLock(void) const;

    //Return status of the pessimistic lock - true for locked, false for unlocked
    bool pessimisticLock(void) const;

    //Return the last in-phase error value
    double lastInPhase(void) const;

    //Return the last quad-phase error value
    double lastQuadPhase(void) const;

    //Return the count on the pessimistic detector
    unsigned pessimisticCount(void) const;

    //Return the count on the optimistic detector
    unsigned optimisticCount(void) const;

private:
    LowPassFilter inPhaseLPF_;

    LowPassFilter quadPhaseLPF_;

    unsigned pessimisticCount_;

    unsigned optimisticCount_;

    bool optimisticLock_;

    bool pessimisticLock_;
};

#endif /* SRC_LOCKDETECTOR_H_ */
