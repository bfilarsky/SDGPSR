#ifndef SRC_NAVBITEDGEDETECTOR_H_
#define SRC_NAVBITEDGEDETECTOR_H_

#include <list>

/*
 * Class for detecting the Nav-Bit edges. The C/A Code repeats every 1 ms, while the Nav-Bits are 20 ms long. When the Nav-Bit is
 * -1, it flips all of the C/A code bits. In order to integrate over the length of a Nav-Bit, we need to know where it starts
 * (it will always be aligned with the C/A sequence). In order to reliably detect the edge, candidates are passed to this class
 * with number of integrations from some reference point. Once a threshold number of candidates are all separated by multiples of 20 ms,
 * the Nav-Bit is determined to have been found
 *
 * | C/A | C/A | C/A | C/A | C/A | C/A | C/A | C/A | C/A | C/A | C/A | C/A | C/A | C/A | C/A | C/A | C/A | C/A | C/A | C/A |
 * |                                                         Nav-Bit                                                       |
 */

class NavBitEdgeDetector {
public:
    NavBitEdgeDetector(unsigned threshold);

    virtual ~NavBitEdgeDetector();

    bool push_back(size_t candidate);

    void clear(void);

private:
    unsigned threshold_;

    std::list<size_t> navBitCandidates_;
};

#endif /* SRC_NAVBITEDGEDETECTOR_H_ */
