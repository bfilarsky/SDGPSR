#ifndef SRC_CACODE_H_
#define SRC_CACODE_H_

/*
 * Coarse/Acquisition Code Generator. On construction, this class will generate the 1023-bit C/A sequence for the selected prn
 * It can then be indexed using the [] operator, and returns the code as +1/-1
 */

#include <iostream>
#include <vector>

const double CHIP_RATE = 1.023e6;
const double CHIP_TIME = 1.0 / CHIP_RATE;
const double CA_CODE_TIME = 1e-3;
const int CA_CODE_LENGTH = 1023;
const unsigned PRN_COUNT = 32;

//Delay taps for PRN-1, as defined in IS-GPS-200
const int G2_DELAY[] = { 5, 6, 7, 8, 17, 18, 139, 140, 141, 251, 252, 254, 255, 256, 257, 258, 469, 470, 471, 472, 473,
        474, 509, 512, 513, 514, 515, 516, 859, 860, 861, 862 };

class CaCode {
public:
    CaCode(unsigned prn);

    virtual ~CaCode();

    int operator[](size_t idx) const;

private:
    unsigned prn_;

    std::vector<int> goldCode_;
};

#endif /* SRC_CACODE_H_ */
