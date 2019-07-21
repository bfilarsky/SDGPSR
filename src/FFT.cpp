#include "FFT.h"

FFT::FFT(unsigned size) {
    fftwVector in(size);
    fftwVector out(size);
    forward_ = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex*>(&in[0]), reinterpret_cast<fftw_complex*>(&out[0]),
            FFTW_FORWARD, FFTW_MEASURE);
    reverse_ = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex*>(&in[0]), reinterpret_cast<fftw_complex*>(&out[0]),
            FFTW_BACKWARD, FFTW_MEASURE);
    forwardIp_ = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex*>(&in[0]),
            reinterpret_cast<fftw_complex*>(&in[0]), FFTW_FORWARD, FFTW_MEASURE);
    reverseIp_ = fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex*>(&in[0]),
            reinterpret_cast<fftw_complex*>(&in[0]), FFTW_BACKWARD, FFTW_MEASURE);
}

FFT::~FFT() {

}

void FFT::forward(std::complex<double> *in, std::complex<double> *out) {
    fftw_execute_dft(forward_, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
}

void FFT::reverse(std::complex<double> *in, std::complex<double> *out) {
    fftw_execute_dft(reverse_, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out));
}

void FFT::forward(std::complex<double> *inplace) {
    fftw_execute_dft(forwardIp_, reinterpret_cast<fftw_complex*>(inplace), reinterpret_cast<fftw_complex*>(inplace));
}

void FFT::reverse(std::complex<double> *inplace) {
    fftw_execute_dft(reverseIp_, reinterpret_cast<fftw_complex*>(inplace), reinterpret_cast<fftw_complex*>(inplace));
}
