/*
 * CodeLockDetector.h
 *
 *  Created on: Apr 8, 2019
 *      Author: bfilarsky
 */

#ifndef SRC_CODELOCKDETECTOR_H_
#define SRC_CODELOCKDETECTOR_H_

#include "LowPassFilter.h"

class CodeLockDetector {
public:
	CodeLockDetector();

	virtual ~CodeLockDetector();

	void error(double magEarly, double magPrompt, double magLate);

	bool optimisticLock(void);

	double early(void){return earlyLPF_.last();}

	double prompt(void){return promptLPF_.last();}

	double late(void){return lateLPF_.last();}

private:
	LowPassFilter earlyLPF_;
	LowPassFilter promptLPF_;
	LowPassFilter lateLPF_;

	unsigned pCount_;
};

#endif /* SRC_CODELOCKDETECTOR_H_ */
