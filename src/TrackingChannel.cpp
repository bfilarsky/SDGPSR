#include "TrackingChannel.h"
#include <unistd.h>

TrackingChannel::TrackingChannel(double fs, unsigned prn, SearchResult searchResult) {
	for (int i = -4; i <= 4; ++i)
		signalTrackers_.push_back(new SignalTracker(fs, prn, searchResult, i*500));
	inputPackets_ = 0;
	prn_ = prn;
}

TrackingChannel::~TrackingChannel() {
	while (signalTrackers_.size()){
		delete signalTrackers_.back();
		signalTrackers_.pop_back();
	}
}

unsigned TrackingChannel::prn(void){
	return prn_;
}

double TrackingChannel::transmitTime(void){
	return signalTrackers_.front()->transmitTime();
}

Vector3d TrackingChannel::satellitePosition(double timeOfWeek){
	return signalTrackers_.front()->satellitePosition(timeOfWeek);
}

complex<double> TrackingChannel::latLong(double timeOfWeek){
	return signalTrackers_.front()->latLong(timeOfWeek);
}

State TrackingChannel::state(void){
	State state = lossOfLock;
	for (auto &tracker : signalTrackers_)
		state = std::max(state, tracker->state());
	return state;
}

bool TrackingChannel::processSamples(fftwVector trackingData){
	for (auto it = signalTrackers_.begin(); it != signalTrackers_.end();){
		if (!(*it)->processSamples(trackingData)){
			cout << "PRN " << prn_ << " channel lost track (" << signalTrackers_.size() - 1 << " left)" << endl;
			delete *it;
			it = signalTrackers_.erase(it);
		}
		else
			++it;
	}

	++inputPackets_;

	if (signalTrackers_.size() > 1 && state() == fullTrack){
		sync();
		for (auto it = signalTrackers_.begin(); it != signalTrackers_.end();){
			if ((*it)->state() != fullTrack){
				delete *it;
				it = signalTrackers_.erase(it);
			}
			else
				++it;
		}
		cout << "PRN " << prn_ << " in full track" << endl;
		while (signalTrackers_.size() > 1){
			delete signalTrackers_.back();
			signalTrackers_.pop_back();
		}
	}
	return signalTrackers_.size();
}

void TrackingChannel::sync(void){
	for (auto &tracker : signalTrackers_)
		tracker->sync();
}
