#ifndef LNAV_WORD_H_
#define LNAV_WORD_H_

#include <vector>
#include <deque>
#include <iostream>

const unsigned WORD_SIZE = 30;
const unsigned SRC_BITS = 24;

/*
 * Class for processing an individual LNAV word. Each word is 30 bits, with 24 bits of data followed by a 6 bit checksum.
 * The checksum is based on the bits of this word, along with bits 29 and 30 (1-indexed) of the previous word. If bit 30
 * of the previous word was 1, all of the bits on this word are inverted.
 */

class LNAV_Word {
public:
    LNAV_Word(const std::deque<int> &word, unsigned lastWordBit29, unsigned lastWordBit30);

    ~LNAV_Word();

    //Word has a valid checksum
    bool valid(void) const;

    //Value of 1-indexed bit
    unsigned bit(size_t idx) const;

    //Interpret bit sequence as unsigned integer. Start is 1-indexed
    unsigned getUnsignedValue(unsigned start, unsigned length) const;

    //Interpret bit sequence as signed integer. Start is 1-indexed
    int getSignedValue(unsigned start, unsigned length) const;

private:
    bool computeChecksum(unsigned lastWordBit29, unsigned lastWordBit30);

    std::vector<unsigned> word_;

    bool valid_;
};

#endif /* LNAV_WORD_H_ */
