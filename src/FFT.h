#ifndef SRC_FFT_H_
#define SRC_FFT_H_

#include <fftw3.h>
#include <vector>
#include <complex>
#include "fftw_allocator.h"

//Standard vector, except that it uses the fftw_allocator, which allocates memory that is aligned for SIMD acceleration
typedef std::vector<std::complex<double>, fftw_allocator<std::complex<double>>> fftwVector;

/*
 * This class is just a container for the fftw3 library. It sets up the library on construction as needed for SDGPSR
 * The functions allow FFTs and IFFTs to be completed both in-place and into a different location
 */

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
