#ifndef SRC_LOWPASSFILTER_H_
#define SRC_LOWPASSFILTER_H_

/*
 * Simple single-pole infinite impulse response low pass filter. General form is y[n] = (1-k)*x[n] + k*y[n-1], where fc = -ln(k)/2pi, or k = exp(-2*pi*fc)
 * Redefining k as 1-k gives y[n] = k*x[n] + (1-k)*y[n-1], k = 1 - exp(-2*pi*fc)
 * Rearranging gives y[n] = k(x[n] - y[n-1]) + y[n-1]
 */

class LowPassFilter {
public:
    LowPassFilter(double fCutoff); //fCutoff must be in range (0,0.5]

    virtual ~LowPassFilter();

    double iterate(double val);    //Add new value and return filtered value

    double last(void) const;       //Last input value

private:
    double k_;

    double last_;//y[n-1]
};

#endif /* SRC_LOWPASSFILTER_H_ */
