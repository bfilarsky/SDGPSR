#ifndef SRC_SNRESTIMATOR_H_
#define SRC_SNRESTIMATOR_H_

#include <complex>

using std::complex;

class SNR_Estimator {
public:
	SNR_Estimator();
	virtual ~SNR_Estimator();

	void input(complex<double> carrier);

	double estimate(void);

private:
	double norm2_;
	double norm4_;
	unsigned sampleCount_;
	double lastEstimate_;
};

#endif /* SRC_CNOESTIMATOR_H_ */
