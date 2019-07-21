#include "OrbitalData.h"
#include <cmath>
#include <iomanip>

const double F_RELATIVISTIC = -4.442807633e-10;

OrbitalData::OrbitalData() {
    currentTow_ = -1;
}

OrbitalData::~OrbitalData() {

}

double OrbitalData::currentNavWordTimeOfWeek(void) {
    if (currentTow_ == -1)
        return -1.0;
    double navBitTime = currentTow_ + WORD_PERIOD * words_.size();
    double timeFromEpochClock = navBitTime - clockData_.Toc_;
    if (timeFromEpochClock > 302400.0)
        timeFromEpochClock -= 604800.0;
    if (timeFromEpochClock < -302400.0)
        timeFromEpochClock += 604800.0;
    double clockCorrection = clockData_.af0_ + clockData_.af1_ * timeFromEpochClock
            + clockData_.af2_ * timeFromEpochClock * timeFromEpochClock;

    double timeFromEpochEphem = navBitTime - currentEphemeris_.Toe_;
    if (timeFromEpochEphem > 302400.0)
        timeFromEpochEphem -= 604800.0;
    if (timeFromEpochEphem < -302400.0)
        timeFromEpochEphem += 604800.0;
    double meanAnomaly = currentEphemeris_.meanAnomaly_ + currentEphemeris_.corrMeanMotion_ * timeFromEpochEphem;
    //Iteratively solve for eccentric anomaly using Newton's method
    double eccentricAnomaly = meanAnomaly;
    eccentricAnomaly -= (eccentricAnomaly - currentEphemeris_.ecc_ * sin(eccentricAnomaly) - meanAnomaly)
            / (1.0 - currentEphemeris_.ecc_ * cos(eccentricAnomaly));
    eccentricAnomaly -= (eccentricAnomaly - currentEphemeris_.ecc_ * sin(eccentricAnomaly) - meanAnomaly)
            / (1.0 - currentEphemeris_.ecc_ * cos(eccentricAnomaly));

    double deltaTr = F_RELATIVISTIC * currentEphemeris_.ecc_ * currentEphemeris_.semiMajorAxisSqrt_
            * sin(eccentricAnomaly);
    double deltaTsv = 0.0; //clockCorrection + deltaTr;
    return navBitTime - deltaTsv;
}

bool OrbitalData::process(const LNAV_Word &word) {
    if (!word.valid()) {
        words_.clear();
        currentTow_ = -1;
        return false;
    }

    words_.push_back(word);

    //All the info we need is on SubFrames 1-3, ignoring SFs 4 & 5
    //See figure 20-2 in IS-GPS-200J
    if (words_.size() == 10) {
        currentTow_ = words_[1].getUnsignedValue(1, 17) * 6;
        switch (words_[1].getUnsignedValue(20, 3)) {
            case 1:
                processSubframe1(words_);
                break;
            case 2:
                processSubframe2(words_);
                break;
            case 3:
                processSubframe3(words_);
                break;
            default:
                break;
        }
        words_.clear();
    }
    return true;
}

//Clock data is all on subframe 1, we can just update all clock info when we get this
//See Figure 20-1 (sheet 1 of 11), and Table 20-I in IS-GPS-200J
void OrbitalData::processSubframe1(const vector<LNAV_Word> &words) {
    clockData_.weekNumber_ = words[2].getUnsignedValue(1, 10);        //Wn
    clockData_.IODC_ = (words[2].getUnsignedValue(23, 2) << 8) | words[6].getUnsignedValue(1, 8); //IODC 8 LSBs
    clockData_.Toc_ = words[7].getUnsignedValue(9, 16) * pow(2, 4);    //Toc
    clockData_.af2_ = words[8].getSignedValue(1, 8) * pow(2.0, -55.0); //Af2
    clockData_.af1_ = words[8].getSignedValue(9, 16) * pow(2.0, -43.0); //Af1
    clockData_.af0_ = words[9].getSignedValue(1, 22) * pow(2.0, -31.0); //Af0
    clockData_.valid_ = true;
}

//Ephemeris is split between two subframes, both have IODE. Set what we can of nextEphemeris here,
//don't update current Ephemeris data until we have everything. We do not want to attempt orbital calculations
//with half our data from one data set (and therefor Time of Validity), and half from another
//See Figure 20-1 (sheet 2 of 11), Table 20-II, and Table 20-III in IS-GPS-200J
void OrbitalData::processSubframe2(const vector<LNAV_Word> &words) {
    nextEphemeris_.IODE_ = words[2].getUnsignedValue(1, 8);
    nextEphemeris_.Crs_ = words[2].getSignedValue(9, 16) * pow(2.0, -5.0);
    nextEphemeris_.meanAnomaly_ = ((words[3].getSignedValue(17, 8) << 24) | (0xffffff & words[4].getSignedValue(1, 24))) * pow(2.0,-31.0) * SEMI_CIRCLE_TO_RAD; /**/
    nextEphemeris_.Cuc_ = words[5].getSignedValue(1, 16) * pow(2.0, -29.0);
    nextEphemeris_.ecc_ = ((words[5].getUnsignedValue(17,8) << 24) | words[6].getUnsignedValue(1,24)) * pow(2.0,-33.0);
    nextEphemeris_.Cus_ = words[7].getSignedValue(1, 16) * pow(2.0, -29.0);
    nextEphemeris_.semiMajorAxisSqrt_ = ((words[7].getUnsignedValue(17,8) << 24) | words[8].getUnsignedValue(1,24)) * pow(2.0,-19.0);
    nextEphemeris_.semiMajorAxis_ = nextEphemeris_.semiMajorAxisSqrt_ * nextEphemeris_.semiMajorAxisSqrt_;
    double n0 = sqrt(MU / (pow(nextEphemeris_.semiMajorAxis_, 3)));
    nextEphemeris_.corrMeanMotion_ = n0 + words[3].getSignedValue(1, 16) * pow(2.0, -43.0) * SEMI_CIRCLE_TO_RAD;
    nextEphemeris_.Toe_ = words[9].getUnsignedValue(1, 16) * pow(2, 4);
}

//This should finish out orbital data we need. Once we have it all, if the IODE in this SF matches the IODE we set in SF2,
//and it is more recent than the current data, update the current data with the new data
//See Figure 20-1 (sheet 3 of 11), Table 20-II, and Table 20-III in IS-GPS-200J
void OrbitalData::processSubframe3(const vector<LNAV_Word> &words) {
    nextEphemeris_.Cic_ = words[2].getSignedValue(1, 16) * pow(2.0, -29.0);
    nextEphemeris_.rightAscension_ = ((words[2].getSignedValue(17, 8) << 24) | (0xffffff & words[3].getSignedValue(1, 24))) * pow(2.0,-31.0) * SEMI_CIRCLE_TO_RAD;
    nextEphemeris_.Cis_ = words[4].getSignedValue(1, 16) * pow(2.0, -29.0);
    nextEphemeris_.inclination_ = ((words[4].getSignedValue(17, 8) << 24) | (0xffffff & words[5].getSignedValue(1, 24))) * pow(2.0,-31.0) * SEMI_CIRCLE_TO_RAD;
    nextEphemeris_.Crc_ = words[6].getSignedValue(1, 16) * pow(2.0, -5.0);
    nextEphemeris_.argOfPerigee_ = ((words[6].getSignedValue(17, 8) << 24) | (0xffffff & words[7].getSignedValue(1, 24))) * pow(2.0,-31.0) * SEMI_CIRCLE_TO_RAD;
    nextEphemeris_.rateOfRightAscen_ = words[8].getSignedValue(1, 24) * pow(2.0, -43.0) * SEMI_CIRCLE_TO_RAD;
    nextEphemeris_.inclinationRate_ = words[9].getSignedValue(9, 14) * pow(2.0, -43.0) * SEMI_CIRCLE_TO_RAD;
    if (nextEphemeris_.IODE_ == words[9].getUnsignedValue(1, 8))
        nextEphemeris_.valid_ = true;
    if (nextEphemeris_.valid_ && (!currentEphemeris_.valid_ || currentEphemeris_.IODE_ != nextEphemeris_.IODE_))
        currentEphemeris_ = nextEphemeris_;
}

//Calculate position of satellite in ECEF at given time of week (sec) //See Table 20-IV in IS-GPS-200J
Vector3d OrbitalData::satellitePosition(double timeOfWeek) {
    Vector3d positionECEF(0.0, 0.0, 0.0);
    if (!currentEphemeris_.valid_)
        return positionECEF;
    double timeFromEpoch = timeOfWeek - currentEphemeris_.Toe_;
    if (timeFromEpoch > 302400.0)
        timeFromEpoch -= 604800.0;
    if (timeFromEpoch < -302400.0)
        timeFromEpoch += 604800.0;

    double meanAnomaly = currentEphemeris_.meanAnomaly_ + currentEphemeris_.corrMeanMotion_ * timeFromEpoch;
    //Iteratively solve for eccentric anomaly using Newton's method
    double eccentricAnomaly = meanAnomaly;
    eccentricAnomaly -= (eccentricAnomaly - currentEphemeris_.ecc_ * sin(eccentricAnomaly) - meanAnomaly)
            / (1.0 - currentEphemeris_.ecc_ * cos(eccentricAnomaly));
    eccentricAnomaly -= (eccentricAnomaly - currentEphemeris_.ecc_ * sin(eccentricAnomaly) - meanAnomaly)
            / (1.0 - currentEphemeris_.ecc_ * cos(eccentricAnomaly));

    double numTrueAnomaly = sqrt(1.0 - pow(currentEphemeris_.ecc_, 2.0)) * sin(eccentricAnomaly);
    double denTrueAnomaly = (cos(eccentricAnomaly) - currentEphemeris_.ecc_);
    double trueAnomaly = atan2(numTrueAnomaly, denTrueAnomaly);

    double argOfLat = trueAnomaly + currentEphemeris_.argOfPerigee_;
    double sin2ArgOfLat = sin(2.0 * argOfLat);
    double cos2ArgOfLat = cos(2.0 * argOfLat);

    double argOfLatCorr = currentEphemeris_.Cus_ * sin2ArgOfLat + currentEphemeris_.Cuc_ * cos2ArgOfLat;
    double radiusCorr = currentEphemeris_.Crs_ * sin2ArgOfLat + currentEphemeris_.Crc_ * cos2ArgOfLat;
    double inclinCorr = currentEphemeris_.Cis_ * sin2ArgOfLat + currentEphemeris_.Cic_ * cos2ArgOfLat;

    double corrArgOfLat = argOfLat + argOfLatCorr;
    double corrRadius = currentEphemeris_.semiMajorAxis_ * (1.0 - currentEphemeris_.ecc_ * cos(eccentricAnomaly))
            + radiusCorr;
    double corrInclin = currentEphemeris_.inclination_ + currentEphemeris_.inclinationRate_ * timeFromEpoch
            + inclinCorr;

    double xkPrime = corrRadius * cos(corrArgOfLat);
    double ykPrime = corrRadius * sin(corrArgOfLat);

    double corrRightAscension = currentEphemeris_.rightAscension_
            + (currentEphemeris_.rateOfRightAscen_ - OMEGA_EARTH) * timeFromEpoch
            - OMEGA_EARTH * currentEphemeris_.Toe_;

    positionECEF.x() = xkPrime * cos(corrRightAscension) - ykPrime * cos(corrInclin) * sin(corrRightAscension);
    positionECEF.y() = xkPrime * sin(corrRightAscension) + ykPrime * cos(corrInclin) * cos(corrRightAscension);
    positionECEF.z() = ykPrime * sin(corrInclin);

    return positionECEF;
}
