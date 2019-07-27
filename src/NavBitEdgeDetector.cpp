#include "NavBitEdgeDetector.h"

NavBitEdgeDetector::NavBitEdgeDetector(unsigned threshold) {
    threshold_ = threshold;
}

NavBitEdgeDetector::~NavBitEdgeDetector() {

}

bool NavBitEdgeDetector::push_back(size_t candidate) {
    navBitCandidates_.push_back(candidate);

    //Once the history is long enough, evaluate the data
    if (navBitCandidates_.size() >= threshold_) {
        //If any of the candidates are not a multiple of 20 integration periods from the first candidate,
        //remove the first candidate and return false. If they all are, clear out the data and return true
        for (auto &candidate : navBitCandidates_) {
            if ((candidate - navBitCandidates_.front()) % 20) {
                navBitCandidates_.pop_front();
                return false;
            }
        }
        navBitCandidates_.clear();
        return true;
    }

    return false;
}

void NavBitEdgeDetector::clear(void) {
    navBitCandidates_.clear();
}
