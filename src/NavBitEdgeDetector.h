/*
 * NavBitEdgeDetector.h
 *
 *  Created on: Apr 7, 2019
 *      Author: bfilarsky
 */

#ifndef SRC_NAVBITEDGEDETECTOR_H_
#define SRC_NAVBITEDGEDETECTOR_H_

#include <deque>

class NavBitEdgeDetector {
public:
	NavBitEdgeDetector();

	virtual ~NavBitEdgeDetector();

	bool push_back(size_t candidate);

	void clear(void);

private:
	std::deque<size_t> navBitCandidates_;
};

#endif /* SRC_NAVBITEDGEDETECTOR_H_ */
