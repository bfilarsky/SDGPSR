#ifndef SRC_FFT_H_
#define SRC_FFT_H_

#include <fftw3.h>
#include <vector>
#include <complex>
#include "fftw_allocator.h"

typedef std::vector<std::complex<double>, fftw_allocator<std::complex<double>>> fftwVector;

class FFT {
public:
    FFT(unsigned size);

    virtual ~FFT();

    void forward(std::complex<double> *in, std::complex<double> *out);

    void reverse(std::complex<double> *in, std::complex<double> *out);

    void forward(std::complex<double> *inplace);

    void reverse(std::complex<double> *inplace);

private:
    fftw_plan forward_;
    fftw_plan reverse_;
    fftw_plan forwardIp_;
    fftw_plan reverseIp_;
};

#endif /* SRC_FFT_H_ */
