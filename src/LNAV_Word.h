/*
 * LNAVWord.h
 *
 *  Created on: Mar 16, 2019
 *      Author: bfilarsky
 */

#ifndef LNAV_WORD_H_
#define LNAV_WORD_H_

#include <vector>
#include <deque>
#include <iostream>

using namespace std;

const unsigned WORD_SIZE = 30;
const unsigned SRC_BITS = 24;

class LNAV_Word {
public:
	LNAV_Word(const deque<int> &word, unsigned lastWordBit29, unsigned lastWordBit30, bool inverted = false);

	~LNAV_Word();

	bool valid(void) const;

	unsigned operator [](size_t idx) const;

	unsigned getUnsignedValue(unsigned start, unsigned length) const;

	int getSignedValue(unsigned start, unsigned length) const;

private:
	bool computeChecksum(unsigned lastWordBit29, unsigned lastWordBit30);

	vector<unsigned> word_;

	bool valid_;
};

#endif /* LNAV_WORD_H_ */
