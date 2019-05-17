#include "LNAV_Word.h"

LNAV_Word::LNAV_Word(const deque<int> &word, unsigned lastWordBit29, unsigned lastWordBit30, bool inverted) {
	word_.resize(WORD_SIZE);
	for (unsigned i = 0; i < WORD_SIZE; ++i){
		word_[i] = (word[i] == 1) ^ inverted;
		if (i < SRC_BITS)
			word_[i] ^= lastWordBit30;
	}
	valid_ = computeChecksum(lastWordBit29, lastWordBit30);
}

LNAV_Word::~LNAV_Word() {

}

bool LNAV_Word::valid(void) const {
	return valid_;
}

unsigned LNAV_Word::operator [](size_t idx) const {
	return word_[idx];
}

unsigned LNAV_Word::getUnsignedValue(unsigned start, unsigned length) const {
	unsigned val = 0;
	for (unsigned i = 0; i < length; ++i)
		val |= word_[i + start - 1] << (length - i - 1);
	return val;
}

int LNAV_Word::getSignedValue(unsigned start, unsigned length) const {
	int val = 0;
	if (word_[start - 1])
	   for (unsigned i = 0; i < 32-length; ++i)
	      val |= 1 << (31 - i);
	for (unsigned i = 0; i < length; ++i){
	   val |= word_[i + start - 1] << (length - i - 1);
	}

	return val;
}

bool LNAV_Word::computeChecksum(unsigned lastWordBit29, unsigned lastWordBit30){
	unsigned bit25 = lastWordBit29 ^ word_[0] ^ word_[1] ^ word_[2] ^ word_[4] ^ word_[5] ^ word_[9] ^ word_[10] ^ word_[11] ^ word_[12] ^
			word_[13] ^ word_[16] ^ word_[17] ^ word_[19] ^ word_[22];
	if (bit25 != word_[24]){
		return false;
	}
	unsigned bit26 = lastWordBit30 ^ word_[1] ^ word_[2] ^ word_[3] ^ word_[5] ^ word_[6] ^ word_[10] ^ word_[11] ^ word_[12] ^ word_[13] ^
			word_[14] ^ word_[17] ^ word_[18] ^ word_[20] ^ word_[23];
	if (bit26 != word_[25]){
		return false;
	}
	unsigned bit27 = lastWordBit29 ^ word_[0] ^ word_[2] ^ word_[3] ^ word_[4] ^ word_[6] ^ word_[7] ^ word_[11] ^ word_[12] ^ word_[13] ^
			word_[14] ^ word_[15] ^ word_[18] ^ word_[19] ^ word_[21];
	if (bit27 != word_[26]){
		return false;
	}
	unsigned bit28 = lastWordBit30 ^ word_[1] ^ word_[3] ^ word_[4] ^ word_[5] ^ word_[7] ^ word_[8] ^ word_[12] ^ word_[13] ^ word_[14] ^
			word_[15] ^ word_[16] ^ word_[19] ^ word_[20] ^ word_[22];
	if (bit28 != word_[27]){
		return false;
	}
	unsigned bit29 = lastWordBit30 ^ word_[0] ^ word_[2] ^ word_[4] ^ word_[5] ^ word_[6] ^ word_[8] ^ word_[9] ^ word_[13] ^ word_[14] ^
			word_[15] ^ word_[16] ^ word_[17] ^ word_[20] ^ word_[21] ^ word_[23];
	if (bit29 != word_[28]){
		return false;
	}
	unsigned bit30 = lastWordBit29 ^ word_[2] ^ word_[4] ^ word_[5] ^ word_[7] ^ word_[8] ^ word_[9] ^ word_[10] ^ word_[12] ^
			word_[14] ^ word_[18] ^ word_[21] ^ word_[22] ^ word_[23];
	if (bit30 != word_[29]){
		return false;
	}
	return true;
}
