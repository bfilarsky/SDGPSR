#include "NavBitEdgeDetector.h"

NavBitEdgeDetector::NavBitEdgeDetector() {

}

NavBitEdgeDetector::~NavBitEdgeDetector() {

}

bool NavBitEdgeDetector::push_back(size_t candidate) {
    navBitCandidates_.push_back(candidate);

    if (navBitCandidates_.size() > 4) {
        for (unsigned i = 1; i < navBitCandidates_.size(); ++i) {
            if ((navBitCandidates_[i] - navBitCandidates_[0]) % 20) {
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
