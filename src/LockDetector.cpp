/*
 * LockDetector.cpp
 *
 *  Created on: Apr 8, 2019
 *      Author: bfilarsky
 */

#include "LockDetector.h"

const double K1 = 0.0247;
const double K2 = 1.5;
const unsigned LP = 50;
const unsigned LO = 240;

LockDetector::LockDetector() : inPhaseLPF_(K1), quadPhaseLPF_(K1) {
	pCount1_ = 0;
	pCount2_ = 0;
	optimisticLock_ = false;
	pessimisticLock_ = false;
}

LockDetector::~LockDetector() {

}

void LockDetector::error(complex<double> error){
	double A = inPhaseLPF_.iterate(fabs(error.real())) / K2;
	double B = quadPhaseLPF_.iterate(fabs(error.imag()));
	if (A > B){
		optimisticLock_ = true;
		pCount2_ = 0;
		if (pCount1_ > LP)
			pessimisticLock_ = true;
		else
			++pCount1_;
	}
	else {
		pessimisticLock_ = false;
		pCount1_ = 0;
		if (pCount2_ > LO)
			optimisticLock_ = false;
		else
			++pCount2_;
	}
}

bool LockDetector::optimisticLock(void){
	return optimisticLock_;
}

bool LockDetector::pessimisticLock(void){
	return pessimisticLock_;
}

double LockDetector::lastInPhase(void){
	return inPhaseLPF_.last();
}

double LockDetector::lastQuadPhase(void){
	return quadPhaseLPF_.last();
}

unsigned LockDetector::pCount1(void){
	return pCount1_;
}

unsigned LockDetector::pCount2(void){
	return pCount2_;
}
