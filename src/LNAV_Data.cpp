#include "LNAV_Data.h"
#include <array>

//The preamble on Word 1 of every subframe. The TLM word is detected by finding this sequence, then verifying the checksum
//on the candidate word
const array<int, 8> TLM_PREAMBLE = { 1, -1, -1, -1, 1, -1, 1, 1 };

LNAV_Data::LNAV_Data() {
    trackingSubFrame_ = false;
    lastWordBit29_ = 0;
    lastWordBit30_ = 0;
}

LNAV_Data::~LNAV_Data() {

}

double LNAV_Data::timeOfLastNavBit(void) const {
    double wordTimeOfWeek = orbitalData_.currentNavWordTimeOfWeek();
    if (wordTimeOfWeek == -1.0)
        return -1.0;
    return wordTimeOfWeek + bitHistory_.size() * BIT_PERIOD;
}

Vector3d LNAV_Data::satellitePosition(double timeOfWeek) const {
    return orbitalData_.satellitePosition(timeOfWeek);
}

bool LNAV_Data::valid(void) const{
    return orbitalData_.valid();
}

void LNAV_Data::navBit(int bit, bool &flipBits) {
    flipBits = false;
    bitHistory_.push_back(bit);
    if (bitHistory_.size() == 30) {
        if (!trackingSubFrame_) {
            if (isTlmWord(bitHistory_, flipBits)) {
                LNAV_Word word(bitHistory_, 0, 0);
                lastWordBit29_ = word.bit(29);
                lastWordBit30_ = word.bit(30);
                bitHistory_.clear();
                trackingSubFrame_ = true;
                orbitalData_.process(word);
            } else
                bitHistory_.pop_front();
        } else {
            LNAV_Word word(bitHistory_, lastWordBit29_, lastWordBit30_);
            lastWordBit29_ = bitHistory_[28] == 1;
            lastWordBit30_ = bitHistory_[29] == 1;
            bitHistory_.clear();
            orbitalData_.process(word);
        }
    }
}

//See Figure 20-2 in IS-GPS-200J
bool LNAV_Data::isTlmWord(std::deque<int> &bits, bool &flipBits) {
    //Check for correlation with preamble
    int correlation = 0;
    for (unsigned i = 0; i < TLM_PREAMBLE.size(); ++i)
        correlation += bits[i] * TLM_PREAMBLE[i];
    if (abs(correlation) < 8)
        return false;
    flipBits = correlation < 0;
    if (flipBits)
        for (auto &bit : bits)
            bit *= -1;

    LNAV_Word tlmCandidate(bits, 0, 0);

    if (tlmCandidate.valid()) {
        return true;
    }
    return false;
}

