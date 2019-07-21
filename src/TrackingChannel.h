#ifndef SRC_TRACKING_CHANNEL_H_
#define SRC_TRACKING_CHANNEL_H_

#include <fstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <thread>
#include <mutex>
#include <list>
#include <memory>
#include "SignalTracker.h"

using std::complex;
using std::cout;
using std::endl;

class TrackingChannel {
public:
    TrackingChannel(double fs, unsigned prn, SearchResult searchResult);

    unsigned prn(void);

    bool processSamples(fftwVector trackingData);

    void sync(void);

    double transmitTime(void);

    Vector3d satellitePosition(double timeOfWeek);

    complex<double> latLong(double timeOfWeek);

    State state(void);

    virtual ~TrackingChannel();

private:
    std::list<std::unique_ptr<SignalTracker>> signalTrackers_;

    unsigned inputPackets_;

    unsigned prn_;
};

#endif /* SRC_GPSSATELLITE_H_ */
