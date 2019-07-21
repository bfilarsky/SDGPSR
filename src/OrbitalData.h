#ifndef ORBITALDATA_H_
#define ORBITALDATA_H_

#include <vector>
#include <eigen3/Eigen/Dense>
#include "LNAV_Word.h"

using Eigen::Vector3d;
using namespace std;

const double MU = 3.986005e14;
const double OMEGA_EARTH = 7.2921151467e-5;
const double GPS_PI = 3.1415926535898;
const double SEMI_CIRCLE_TO_RAD = GPS_PI;
const double WORD_PERIOD = 0.6;

struct ClockData {
    unsigned weekNumber_;
    unsigned IODC_;
    unsigned Toc_;
    double Tgd_;
    double af2_;
    double af1_;
    double af0_;
    bool valid_;

    ClockData() {
        weekNumber_ = 0;
        IODC_ = 0;
        Toc_ = 0;
        Tgd_ = 0.0;
        af2_ = 0.0;
        af1_ = 0.0;
        af0_ = 0.0;
        valid_ = false;
    }
};

struct EphemData {
    EphemData() {
        IODE_ = 0;
        Crs_ = 0.0;
        corrMeanMotion_ = 0.0;
        meanAnomaly_ = 0.0;
        Cuc_ = 0.0;
        ecc_ = 0.0;
        Cus_ = 0.0;
        semiMajorAxisSqrt_ = 0.0;
        semiMajorAxis_ = 0.0;
        Toe_ = 0;
        Cic_ = 0.0;
        rightAscension_ = 0.0;
        Cis_ = 0.0;
        inclination_ = 0.0;
        Crc_ = 0.0;
        argOfPerigee_ = 0.0;
        rateOfRightAscen_ = 0.0;
        inclinationRate_ = 0;
        valid_ = false;
    }

    unsigned IODE_;
    double Crs_;
    double corrMeanMotion_;
    double meanAnomaly_;
    double Cuc_;
    double ecc_;
    double Cus_;
    double semiMajorAxis_;
    double semiMajorAxisSqrt_;
    unsigned Toe_;
    double Cic_;
    double rightAscension_;
    double Cis_;
    double inclination_;
    double Crc_;
    double argOfPerigee_;
    double rateOfRightAscen_;
    double inclinationRate_;
    bool valid_;
};

class OrbitalData {
public:
    OrbitalData();
    virtual ~OrbitalData();

    double currentNavWordTimeOfWeek(void);

    bool process(const LNAV_Word &word);

    Vector3d satellitePosition(double timeOfWeek);

private:
    void processSubframe1(const vector<LNAV_Word> &words);

    void processSubframe2(const vector<LNAV_Word> &words);

    void processSubframe3(const vector<LNAV_Word> &words);

    vector<LNAV_Word> words_;

    int currentTow_;

    EphemData currentEphemeris_;

    EphemData nextEphemeris_;

    ClockData clockData_;
};

#endif /* ORBITALDATA_H_ */
