#ifndef SRC_SNRESTIMATOR_H_
#define SRC_SNRESTIMATOR_H_

#include <complex>

/*
 * Estimate the Signal/Noise ratio given an input carrier
 */

class SNR_Estimator {
public:
    SNR_Estimator();

    virtual ~SNR_Estimator();

    void input(std::complex<double> carrier);

    double estimate(void);

private:
    double norm2_;
    double norm4_;
    unsigned sampleCount_;
    double lastEstimate_;
};

#endif /* SRC_CNOESTIMATOR_H_ */
