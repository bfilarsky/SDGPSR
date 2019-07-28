#ifndef SRC_SDGPSR_H_
#define SRC_SDGPSR_H_

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <queue>
#include <vector>
#include <cstring>
#include <thread>
#include <fftw3.h>
#include <list>
#include <memory>
#include <unordered_map>
#include <eigen3/Eigen/Dense>
#include <atomic>
#include "CaCode.h"
#include "FFT.h"
#include "SignalTracker.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector4d;

const uint64_t GPS_L1_HZ = 1575420000;
const double SAT_FOUND_THRESH = 10.0;
const double SPEED_OF_LIGHT_MPS = 299792458.0;

/*
 * Software-Defined Global Positioning System Receiver (SDGPSR) runs at sample rate fs, and takes in baseband IQ data in 1ms intervals
 * using the basebandSignal function. The data should be roughly centered around clockOffset. A minimum of about 35 seconds worth of data
 * is required to get a navigation solution, as each frame takes 30 seconds, and some data is used in the search and tracker initialization
 * processes. SDGPSR launches its own threads to handle calculation, so the functions here return immediately.
 */

class SDGPSR {
public:
    //fs = data clock rate, clock offset is used to inform system of a hardware clock error
    SDGPSR(double fs, double clockOffsetHz);

    virtual ~SDGPSR();

    //Enque data for processing. Processing is handled by threads of this
    //class, so this returns after data is copied into the queue
    void basebandSignal(const fftwVector &data);

    void basebandSignal(fftwVector &&data);

    //Check to see if all of the enqued data has been processed
    bool synced(void);

    //Get user position in WGS84 ECEF frame (m)
    Vector3d positionECEF(void);

    //Get user position in Lat/Lon/Alt (deg,deg,m)
    Vector3d positionLLA(void);

    //Get user GPS time of week (s)
    double timeOfWeek(void);

    std::vector<std::pair<unsigned, State>> trackingStatus();

    bool navSolution(void);

private:
    void solve(void);

    void signalProcessing();

    std::vector<double> nonCoherentCorrelator(std::vector<fftwVector> &searchData, fftwVector &basebandCode, unsigned corrCount);

    void basebandGenerator(unsigned prn, fftwVector &basebandCode, double freqOffset);

    SearchResult search(std::vector<fftwVector> &searchData,
                        unsigned prn,
                        unsigned corrCount,
                        double freqStart,
                        double freqStop,
                        double freqStep);

    Vector4d userEstimateEcefTime_;
    //Mutex should be locked any time a public function reads from userEstimateEcefTime_, and any time it is written to
    std::mutex userEstMutex_;

    bool navSolutionStarted_;

    double fs_;

    std::queue<fftwVector> input_;
    //Mutex should be locked any time a public function reads from input_, and any time it is modified
    std::mutex inputMutex_;

    FFT fft_;

    atomic_bool run_;

    std::list<std::unique_ptr<SignalTracker>> channels_;
    //Mutex should be locked any time a public function reads from channels_, and any time it is modified
    std::mutex channelMutex_;

    double clockOffset_;

    std::thread signalProcessor_;

#ifdef DEBUG_FILES
    ofstream userEstimates_;
    std::unordered_map<unsigned,ofstream> residualsOutput_;
#endif
};

#endif /* SRC_SDGPSR_H_ */
