#ifndef SRC_LNAV_DATA_H_
#define SRC_LNAV_DATA_H_

#include <fstream>
#include <iostream>
#include <iomanip>
#include <deque>
#include <vector>
#include "OrbitalData.h"

const double BIT_PERIOD = 0.02;

using namespace std;

class LNAV_Data {
public:
    LNAV_Data();

    virtual ~LNAV_Data();

    void navBit(int bit, bool &flipBits);

    double timeOfLastNavBit(void);

    Vector3d satellitePosition(double timeOfWeek);

private:
    bool isTlmWord(deque<int> &bits, bool &flipBits);

    deque<int> bitHistory_;

    OrbitalData orbitalData_;

    bool trackingSubFrame_;

    unsigned lastWordBit29_;

    unsigned lastWordBit30_;
};

#endif /* SRC_LNAV_DATA_H_ */
