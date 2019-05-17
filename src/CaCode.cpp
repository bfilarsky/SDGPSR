/*
 * CaCode.cpp
 *
 *  Created on: Feb 9, 2019
 *      Author: bfilarsky
 */

#include "CaCode.h"

CaCode::CaCode(unsigned prn) {
	prn_ = prn;
	std::vector<int> G1(CA_CODE_LENGTH, 1);
	std::vector<int> G2(CA_CODE_LENGTH, 1);
	for (int i = 10; i < CA_CODE_LENGTH; ++i){
		G1[i] = (G1[i-3] + G1[i-10]) % 2;
		G2[i] = (G2[i-2] + G2[i-3] + G2[i-6] + G2[i-8] + G2[i-9] + G2[i-10]) % 2;
	}

	goldCode_.resize(CA_CODE_LENGTH);
	for (int i = 0; i < CA_CODE_LENGTH; ++i){
		int g2Idx = (i - G2_DELAY[prn_-1] + CA_CODE_LENGTH) % CA_CODE_LENGTH;
		goldCode_[i] = (G1[i] + G2[g2Idx]) % 2;
	}
}

CaCode::~CaCode() {

}

int CaCode::operator[](size_t idx) const {
	return goldCode_[idx];
}
