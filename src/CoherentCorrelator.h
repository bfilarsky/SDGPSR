#ifndef SRC_COHERENTCORRELATOR_H_
#define SRC_COHERENTCORRELATOR_H_

#include <iostream>
#include <complex>
#include <fstream>
#include "FFT.h"
#include "CaCode.h"

using std::cout;
using std::endl;

class CoherentCorrelator {
public:
	CoherentCorrelator(double fs, unsigned prn, bool recordOn = false);

	virtual ~CoherentCorrelator();

	double integrate(fftwVector &data, double codeFreqOffset, double timeOffset, unsigned integrationCount, bool *ready = nullptr);

	std::complex<double> dump(void);

	std::complex<double> val(void);

	double integrationTime(void);

private:
	CaCode caCode_;

	double fs_;

	std::complex<double> currentIntegration_;

	std::complex<double> readyVal_;

	bool ready_;

	unsigned integrationCount_;

	int lastCaChipNum_;

	bool recorderOn_;
	std::ofstream recorder_;
};

#endif /* SRC_COHERENTCORRELATOR_H_ */
