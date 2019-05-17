/*
 * LowPassFilter.cpp
 *
 *  Created on: Apr 7, 2019
 *      Author: bfilarsky
 */

#include "LowPassFilter.h"

LowPassFilter::LowPassFilter(double k) {
	k_   = k;
	last_ = 0.0;
}

LowPassFilter::~LowPassFilter() {

}

double LowPassFilter::iterate(double input){
	double retVal = (input - last_) * k_ + last_;
	last_ = retVal;
	return retVal;
}

double LowPassFilter::last(void){
	return last_;
}
