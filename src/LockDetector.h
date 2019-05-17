#ifndef SRC_LOCKDETECTOR_H_
#define SRC_LOCKDETECTOR_H_

#include <complex>
#include "LowPassFilter.h"

using std::complex;

class LockDetector {
public:
	LockDetector();

	virtual ~LockDetector();

	void error(complex<double> error);

	bool optimisticLock(void);

	bool pessimisticLock(void);

	double lastInPhase(void);

	double lastQuadPhase(void);

	unsigned pCount1(void);

	unsigned pCount2(void);

private:
	LowPassFilter inPhaseLPF_;

	LowPassFilter quadPhaseLPF_;

	unsigned pCount1_;

	unsigned pCount2_;

	bool optimisticLock_;

	bool pessimisticLock_;
};

#endif /* SRC_LOCKDETECTOR_H_ */
