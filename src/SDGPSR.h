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
#include <eigen3/Eigen/Dense>
#include <atomic>
#include "CaCode.h"
#include "FFT.h"
#include "TrackingChannel.h"

using std::cout;
using std::endl;
using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector4d;

const uint64_t GPS_L1_HZ = 1575420000;
const double SAT_FOUND_THRESH = 10.0;
const double SPEED_OF_LIGHT_MPS = 299792458;

class SDGPSR {
public:
    SDGPSR(double fs, double clockOffset);

    virtual ~SDGPSR();

    void basebandSignal(fftwVector &data);

    void sync(void);

private:
    void solve(void);

    void signalProcessing();

    std::vector<double> nonCoherentCorrelator(std::vector<fftwVector> &searchData, fftwVector &basebandCode,
            unsigned corrCount);

    void basebandGenerator(unsigned prn, fftwVector &basebandCode, double freqOffset);

    SearchResult search(std::vector<fftwVector> &searchData, unsigned prn, unsigned corrCount, double freqStart,
            double freqStop, double freqStep);

    Vector4d userEstimateEcefTime_;

    double fs_;

    std::queue<fftwVector> input_;

    std::mutex inputMutex_;

    FFT fft_;

    atomic_bool run_;

    std::list<std::unique_ptr<TrackingChannel>> channels_;

    double clockOffset_;

    std::thread signalProcessor_;

    ofstream userEstimates_;
    ofstream innovations_;
};

#endif /* SRC_SDGPSR_H_ */
