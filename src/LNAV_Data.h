#ifndef SRC_LNAV_DATA_H_
#define SRC_LNAV_DATA_H_

#include <deque>
#include <eigen3/Eigen/Dense>
#include "OrbitalData.h"

/*
 * Class for reading in the LNAV bitstream. This class takes bits in, determines the start of a subframe, and passes valid words
 * to the Orbital Data class for processing
 */

const double BIT_PERIOD = 0.02;

class LNAV_Data {
public:
    LNAV_Data();

    virtual ~LNAV_Data();

    //Pass a nav bit in for processing. Flip bits will be set to true if it determines the bits are inverted
    void navBit(int bit, bool &flipBits);

    //GPS Time of Week at the end of the last nav bit
    double timeOfLastNavBit(void) const;

    //ECEF position of the satellite at the given time of week
    Vector3d satellitePosition(double timeOfWeek) const;

    bool valid(void) const;

private:
    bool isTlmWord(std::deque<int> &bits, bool &flipBits);

    std::deque<int> bitHistory_;

    OrbitalData orbitalData_;

    bool trackingSubFrame_;

    unsigned lastWordBit29_;

    unsigned lastWordBit30_;
};

#endif /* SRC_LNAV_DATA_H_ */
