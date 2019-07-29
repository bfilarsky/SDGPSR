#ifndef SRC_COHERENTCORRELATOR_H_
#define SRC_COHERENTCORRELATOR_H_

#include <complex>
#include "FFT.h"
#include "CaCode.h"

/*
 * This is the correlator used by the signal tracker. This computes the complex correlation between the received signal (which should be wiped of carrier & Doppler),
 * and the C/A Code. It correlates in 1ms intervals, and the integration count gives the number of 1ms intervals to run before calling ready.
 */

class CoherentCorrelator {
public:
    CoherentCorrelator(double fs, unsigned prn);

    virtual ~CoherentCorrelator();

    //Add the correlation of data to the currentIntegration value. TimeOffset is the time within the C/A code period (1 ms) to start the correlation
    //If a non-null pointer is passed to ready, it will update the bool with a true when integrationCount has been met
    double integrate(fftwVector &data, double codeFreqOffset, double timeOffset, unsigned integrationCount,
            bool *ready = nullptr);

    //Get the next integrated value and dump it.
    std::complex<double> dump(void);

    //Get the current integrated value, but do not dump it
    std::complex<double> val(void) const;

    //Get the total time the integrator has run
    double integrationTime(void) const;

private:
    CaCode caCode_;

    double fs_;

    std::complex<double> currentIntegration_;

    std::complex<double> readyVal_;

    bool ready_;

    unsigned integrationCount_;

    int lastCaChipNum_;
};

#endif /* SRC_COHERENTCORRELATOR_H_ */
