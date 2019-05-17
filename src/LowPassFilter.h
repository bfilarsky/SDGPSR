/*
 * LowPassFilter.h
 *
 *  Created on: Apr 7, 2019
 *      Author: bfilarsky
 */

#ifndef SRC_LOWPASSFILTER_H_
#define SRC_LOWPASSFILTER_H_

#include <vector>

class LowPassFilter {
public:
	LowPassFilter(double k);

	virtual ~LowPassFilter();

	double iterate(double val);

	double last(void);

private:
	double k_;

	double last_;
};

#endif /* SRC_LOWPASSFILTER_H_ */
