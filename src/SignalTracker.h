#ifndef SRC_SIGNAL_TRACKER_H_
#define SRC_SIGNAL_TRACKER_H_

#include <fstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <thread>
#include <mutex>
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

using std::cout;
using std::endl;

const double CORR_OFFSET = CHIP_TIME/2.0;

enum State {
	lossOfLock,
	closingCarrierFLL,
	closingCarrierPLL,
	findingNavBitEdge,
	fullTrack
};

struct SearchResult{
	bool found;
	double baseBandFreq;
	double power;
	int sampleOffset;

	SearchResult(){
		found = false;
		baseBandFreq = 0.0;
		power        = 0.0;
		sampleOffset = 0;
	}
};

class SignalTracker {
public:
	SignalTracker(double fs, unsigned prn, SearchResult searchResult, double searchFreqOffset);

	unsigned prn(void);

	bool processSamples(fftwVector trackingData);

	void sync(void);

	double transmitTime(void);

	Vector3d satellitePosition(double timeOfWeek);

	complex<double> latLong(double timeOfWeek);

	State state(void);

	double CNoEst(void);

	virtual ~SignalTracker();

private:
	std::complex<double> coherentCorrelator(std::complex<double> *data, double chipFreqOffset, double timeOffset, unsigned intervals);

	void frequencyShift(fftwVector &data, double frequencyShift_Hz);

	void threadFunction(void);

	bool runThread_;

	std::mutex trackingDataAccess_;

	std::list<fftwVector> trackingData_;

	LNAV_Data lnav_data_;

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
	LockDetector carrierFllLock_;
	LockDetector carrierPllLock_;
	LowPassFilter carrierPhaseLPF_;
	LowPassFilter carrierFreqLPF_;

	CodeLockDetector codeLock_;

	unsigned prn_;

	double fs_;

	State state_;

	std::ofstream codeRecorder_;
	std::ofstream carrierRecorder_;

	std::thread processingThread_;
};

#endif /* SRC_GPSSATELLITE_H_ */
