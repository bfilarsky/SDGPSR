#include "CodeLockDetector.h"

const double K1 = 0.0247;
const unsigned LO = 240;

CodeLockDetector::CodeLockDetector() : earlyLPF_(K1), promptLPF_(K1), lateLPF_(K1) {
	pCount_ = 0;
}

CodeLockDetector::~CodeLockDetector() {

}

void CodeLockDetector::error(double magEarly, double magPrompt, double magLate){
	double early  = earlyLPF_.iterate(magEarly);
	double prompt = promptLPF_.iterate(magPrompt);
	double late   = lateLPF_.iterate(magLate);

	if (early > prompt || late > prompt)
		++pCount_;
	else
		pCount_ = 0;
}

bool CodeLockDetector::optimisticLock(void){
	return pCount_ < LO;
}
