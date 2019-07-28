#ifndef SRC_SIGNAL_TRACKER_H_
#define SRC_SIGNAL_TRACKER_H_

#include <fstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <thread>
#include <mutex>
#include <atomic>
#include <list>
#include "CaCode.h"
#include "FFT.h"
#include "CoherentCorrelator.h"
#include "LNAV_Data.h"
#include "SNR_Estimator.h"
#include "LowPassFilter.h"
#include "NavBitEdgeDetector.h"
#include "LockDetector.h"
#include "CodeLockDetector.h"

const double CORR_OFFSET = CHIP_TIME / 2.0;

enum State {
    lossOfLock,
    closingCarrierFLL,
    closingCarrierPLL,
    findingNavBitEdge,
    fullTrack,
    fullNav
};

struct SearchResult {
    bool found;
    double baseBandFreq;
    double power;
    int sampleOffset;

    SearchResult() {
        found = false;
        baseBandFreq = 0.0;
        power = 0.0;
        sampleOffset = 0;
    }
};

/*
 * Signal tracker is the heart of the receiver - one of these tracks each satellite. This class is responsible for the
 * majority of the signal processing work that SDGPSR does. This class maintains its own thread to process all of the
 * data passed to it
 */

class SignalTracker {
public:
    SignalTracker(double fs, unsigned prn, SearchResult searchResult);

    virtual ~SignalTracker();

    //PRN of satellite being tracked
    unsigned prn(void);

    //Copy data into queue for processing
    State processSamples(fftwVector trackingData);

    //Blocking call that returns when the channel has processed all data passed to it
    void sync(void);

    //Transmit time of last sample received
    double transmitTime(void);

    //WGS84 position of the satellite at the given GPS time of week
    Vector3d satellitePosition(double timeOfWeek);

    //Navigation state of signal tracker
    State state(void);

    //Carrier/Noise ratio estimate (linear)
    double CNoEst(void);

private:
    std::complex<double> coherentCorrelator(std::complex<double> *data,
                                            double chipFreqOffset,
                                            double timeOffset,
                                            unsigned intervals);

    void frequencyShift(fftwVector &data, double frequencyShift_Hz);

    void threadFunction(void);

    atomic_bool runThread_;

    atomic_bool synced_;

    std::mutex trackingDataAccess_;

    std::list<fftwVector> trackingData_;

    LNAV_Data lnav_data_;
    std::mutex lnavMutex_;

    CoherentCorrelator early_;
    CoherentCorrelator prompt_;
    CoherentCorrelator late_;

    CoherentCorrelator carrierCorrelator_;
    complex<double> lastCarrier_;

    double timeSinceStart_;
    size_t processedPackets_;

    double carrierPhase_;
    double carrierFreq_;
    double codetime_;
    double codeFreq_;

    unsigned integrationLength_; //(Number of 1 ms intervals to integrate)

    NavBitEdgeDetector navBitEdgeDetector_;

    SNR_Estimator snrEstimator_;
    std::mutex snrMutex_;

    LockDetector carrierFllLock_;
    LockDetector carrierPllLock_;
    LowPassFilter carrierPhaseLPF_;
    LowPassFilter carrierFreqLPF_;

    CodeLockDetector codeLock_;

    unsigned prn_;

    double fs_;

    atomic<State> state_;

#ifdef DEBUG_FILES
    std::ofstream codeRecorder_;
    std::ofstream carrierRecorder_;
#endif

    std::thread processingThread_;
};

#endif /* SRC_GPSSATELLITE_H_ */
