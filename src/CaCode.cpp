#include "CaCode.h"
#include <stdexcept>

CaCode::CaCode(unsigned prn) {
    if (prn == 0 || prn > PRN_COUNT)
        throw std::range_error("PRN out of range");
    prn_ = prn;

    //Compute G1 and G2 per IS-GPS-200
    std::vector<int> G1(CA_CODE_LENGTH, 1);
    std::vector<int> G2(CA_CODE_LENGTH, 1);
    for (int i = 10; i < CA_CODE_LENGTH; ++i) {
        G1[i] = (G1[i - 3] + G1[i - 10]) % 2;
        G2[i] = (G2[i - 2] + G2[i - 3] + G2[i - 6] + G2[i - 8] + G2[i - 9] + G2[i - 10]) % 2;
    }

    //Compute C/A code for selected PRN from G1 and G2 per IS-GPS-200, and set 0s to -1
    goldCode_.resize(CA_CODE_LENGTH);
    for (int i = 0; i < CA_CODE_LENGTH; ++i) {
        int g2Idx = (i - G2_DELAY[prn_ - 1] + CA_CODE_LENGTH) % CA_CODE_LENGTH;
        goldCode_[i] = (G1[i] + G2[g2Idx]) % 2;
        if (goldCode_[i] == 0)
            goldCode_[i] = -1;
    }
}

CaCode::~CaCode() {

}

int CaCode::operator[](size_t idx) const {
    return goldCode_[idx];
}
