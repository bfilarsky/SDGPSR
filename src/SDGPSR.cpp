#include "SDGPSR.h"
#include <iomanip>
#include <unistd.h>

const double SEARCH_WINDOW_BANDWIDTH = 10e3;
const double SEARCH_WINDOW_STEP_SIZE = 500.0;

SDGPSR::SDGPSR(double fs, double clockOffset) : fs_(fs), fft_(fs_ * CA_CODE_TIME), run_(true), signalProcessor_(
        &SDGPSR::signalProcessing, this) {
    clockOffset_ = clockOffset;
    userEstimateEcefTime_ = Vector4d(6378137.0, 0.0, 0.0, 0.0);

    userEstimates_.open("userEstimates.bin", std::ofstream::binary);
    innovations_.open("innovations.bin", std::ofstream::binary);
}

SDGPSR::~SDGPSR() {
    run_ = false;
    signalProcessor_.join();
}

void SDGPSR::basebandSignal(fftwVector &data) {
    inputMutex_.lock();
    input_.push(data);
    inputMutex_.unlock();
}

void SDGPSR::sync(void) {
    while (1) {
        inputMutex_.lock();
        size_t inputSize = input_.size();
        inputMutex_.unlock();
        if (!inputSize)
            return;
        usleep(1e3);
    }
}

std::vector<double> SDGPSR::nonCoherentCorrelator(std::vector<fftwVector> &searchData, fftwVector &basebandCode,
        unsigned corrCount) {
    std::vector<double> result(basebandCode.size());
    fftwVector result_f(basebandCode.size());

    fft_.forward(&basebandCode[0]);

    for (unsigned corrNumber = 0; corrNumber < corrCount; ++corrNumber) {
        for (unsigned sample = 0; sample < basebandCode.size(); ++sample) {
            result_f[sample] = basebandCode[sample] * std::conj(searchData[corrNumber][sample]);
        }
        fft_.reverse(&result_f[0]);
        for (unsigned sample = 0; sample < basebandCode.size(); ++sample) {
            result[sample] += abs(result_f[sample]) / corrCount;
        }
    }

    return result;
}

void SDGPSR::basebandGenerator(unsigned prn, fftwVector &basebandCode, double freqOffset) {
    CaCode caCode(prn);
    basebandCode.resize(fs_ * CA_CODE_TIME);
    for (unsigned i = 0; i < basebandCode.size(); ++i) {
        double time = i / fs_;
        int chip = caCode[time * CHIP_RATE] == 1 ? 1 : -1;
        double phase = 2.0 * M_PI * freqOffset * time;
        basebandCode[i] = std::complex<double>(chip * cos(phase), chip * sin(phase));
    }
}

SearchResult SDGPSR::search(std::vector<fftwVector> &searchData, unsigned prn, unsigned corrCount, double freqStart,
        double freqStop, double freqStep) {
    unsigned freqWindowSize = abs((freqStart - freqStop) / freqStep) + 1;
    std::vector<std::vector<double>> searchWindow(freqWindowSize);
    fftwVector basebandCode;

    double max = 0.0;
    double sum = 0.0;
    double peakFreq = 0.0;
    unsigned sampleOffset = 0;
    for (unsigned freqWindow = 0; freqWindow < freqWindowSize; ++freqWindow) {
        double freq = freqStart + freqWindow * freqStep;
        basebandGenerator(prn, basebandCode, freq);
        searchWindow[freqWindow] = nonCoherentCorrelator(searchData, basebandCode, corrCount);

        for (unsigned i = 0; i < searchWindow[freqWindow].size(); ++i) {
            sum += searchWindow[freqWindow][i];
            if (searchWindow[freqWindow][i] > max) {
                max = searchWindow[freqWindow][i];
                peakFreq = freq;
                sampleOffset = i;
            }
        }
    }

    std::string filename = "prn" + std::to_string(prn) + "searchWindow" + std::to_string((int) freqStep) + ".bin";
    std::ofstream output(filename, std::ofstream::binary);
    for (unsigned i = 0; i < searchWindow.size(); ++i)
        output.write((char*) &searchWindow[i][0], searchWindow[i].size() * sizeof(searchWindow[i][0]));

    double mean = sum / (freqWindowSize * basebandCode.size());
    double accum = 0.0;
    for (unsigned freqWindow = 0; freqWindow < freqWindowSize; ++freqWindow)
        for (unsigned i = 0; i < searchWindow[freqWindow].size(); ++i)
            accum += (searchWindow[freqWindow][i] - mean) * (searchWindow[freqWindow][i] - mean);
    double stdDev = sqrt(accum / (freqWindowSize * basebandCode.size() - 1));
    double maxAboveMeanStdDevs = (max - mean) / stdDev;

    SearchResult searchResult;
    if (maxAboveMeanStdDevs > SAT_FOUND_THRESH) {
        searchResult.found = true;
        searchResult.baseBandFreq = peakFreq;
        searchResult.power = maxAboveMeanStdDevs;
        searchResult.sampleOffset = sampleOffset;
    }

    return searchResult;
}

void SDGPSR::solve(void) {
    vector<Vector4d> satPosAndTime;
    double maxTime = 0.0;
    cout << channels_.front()->transmitTime() << endl;
    for (auto &chan : channels_) {
        chan->sync();
        double time = chan->transmitTime();
        Vector3d position = chan->satellitePosition(time);
        if (userEstimateEcefTime_[3] != 0.0) {
            double tof = userEstimateEcefTime_[3] - time;
            double earthRotation = -tof * OMEGA_EARTH;
            position[0] = position[0] * cos(earthRotation) - position[1] * sin(earthRotation);
            position[1] = position[0] * sin(earthRotation) + position[1] * cos(earthRotation);
        }
        if (position != Vector3d(0.0, 0.0, 0.0)) {
            Vector4d temp;
            temp[0] = position.x();
            temp[1] = position.y();
            temp[2] = position.z();
            temp[3] = time;
            if (time > maxTime)
                maxTime = time;
            satPosAndTime.push_back(temp);
        }
    }

    if (satPosAndTime.size()) {
        if (userEstimateEcefTime_[3] == 0.0)
            userEstimateEcefTime_[3] = maxTime + .08;
        MatrixXd hMatrix(satPosAndTime.size(), 4);
        VectorXd deltaPseudoranges(satPosAndTime.size());
        for (unsigned i = 0; i < satPosAndTime.size(); ++i) {
            double xRangeEst = -(satPosAndTime[i][0] - userEstimateEcefTime_[0]);
            double yRangeEst = -(satPosAndTime[i][1] - userEstimateEcefTime_[1]);
            double zRangeEst = -(satPosAndTime[i][2] - userEstimateEcefTime_[2]);
            double rangeEstimate = sqrt(xRangeEst * xRangeEst + yRangeEst * yRangeEst + zRangeEst * zRangeEst);
            double pseudorangeEst = (userEstimateEcefTime_[3] - satPosAndTime[i][3]) * SPEED_OF_LIGHT_MPS;
            deltaPseudoranges(i) = rangeEstimate - pseudorangeEst;

            double di = i;
            innovations_.write((char*) &di, sizeof(di));
            innovations_.write((char*) &rangeEstimate, sizeof(rangeEstimate));
            innovations_.write((char*) &pseudorangeEst, sizeof(pseudorangeEst));

            hMatrix(i, 0) = xRangeEst / rangeEstimate;
            hMatrix(i, 1) = yRangeEst / rangeEstimate;
            hMatrix(i, 2) = zRangeEst / rangeEstimate;
            hMatrix(i, 3) = 1.0;
        }
        Vector4d deltaEst = (hMatrix.transpose() * hMatrix).ldlt().solve(hMatrix.transpose() * deltaPseudoranges);
        userEstimateEcefTime_[0] -= deltaEst[0];
        userEstimateEcefTime_[1] -= deltaEst[1];
        userEstimateEcefTime_[2] -= deltaEst[2];
        userEstimateEcefTime_[3] += deltaEst[3] / SPEED_OF_LIGHT_MPS;

        userEstimates_.write((char*) &userEstimateEcefTime_, sizeof(userEstimateEcefTime_));

        //Poor man's Lat/Lon conversion - assuming earth is perfect sphere
        double lat = atan2(userEstimateEcefTime_.z(), sqrt(userEstimateEcefTime_.x() * userEstimateEcefTime_.x() + userEstimateEcefTime_.y() * userEstimateEcefTime_.y())) * 180.0 / M_PI;
        double lon = atan2(userEstimateEcefTime_.y(), userEstimateEcefTime_.x()) * 180.0 / M_PI;
        cout << floor(fabs(lat)) << ' ' << floor(fmod(fabs(lat), 1.0) * 60.0) << "\' " << floor(fmod(fmod(fabs(lat), 1.0) * 60.0, 1.0) * 60.0) << "\""<< endl;
        cout << floor(fabs(lon)) << ' ' << floor(fmod(fabs(lon), 1.0) * 60.0) << "\' " << floor(fmod(fmod(fabs(lon), 1.0) * 60.0, 1.0) * 60.0) << "\""<< endl;
        cout << userEstimateEcefTime_ - Vector4d(6378137.0, 0.0, 0.0, 0.0) << endl << endl;
    }
}

void SDGPSR::signalProcessing() {
    //FFT the first CORR_COUNT packets in order to conduct search
    const unsigned CORR_COUNT = 128;
    std::vector<fftwVector> searchData(CORR_COUNT);
    for (unsigned i = 0; i < CORR_COUNT; ++i) {
        while (1) {
            inputMutex_.lock();
            size_t inputSize = input_.size();
            inputMutex_.unlock();
            if (inputSize)
                break;
            if (!run_)
                return;
            usleep(1e3);
        }
        inputMutex_.lock();
        searchData[i].resize(input_.front().size());
        fft_.forward(&input_.front()[0], &searchData[i][0]);
        input_.pop();
        inputMutex_.unlock();
    }

    //Using the FFT'd data, search for all PRNs
    for (unsigned prn = 1; prn <= 32; ++prn) {
        cout << prn << endl;
        SearchResult searchResult = search(searchData, prn, CORR_COUNT, clockOffset_ - SEARCH_WINDOW_BANDWIDTH / 2.0,
                clockOffset_ + SEARCH_WINDOW_BANDWIDTH / 2.0, SEARCH_WINDOW_STEP_SIZE);
        if (searchResult.found) {
            cout << prn << ',' << searchResult.power << ',' << searchResult.sampleOffset << ','
                    << searchResult.baseBandFreq << endl;
            channels_.push_back(std::unique_ptr<TrackingChannel>(new TrackingChannel(fs_, prn, searchResult)));
        }
    }

    //Ensure the minimum number of required satellites have been found
    const unsigned MIN_SATELLITES = 4;
    if (channels_.size() < MIN_SATELLITES) {
        cout << channels_.size() << " satellites found. Minimum " << MIN_SATELLITES << " required" << endl;
        exit(1);
    }

    //Begin tracking satellites
    unsigned trackingPacketCount = 0;
    while (run_) {
        for (auto chanIterator = channels_.begin(); chanIterator != channels_.end();) {
            while (1) {
                inputMutex_.lock();
                size_t inputSize = input_.size();
                inputMutex_.unlock();
                if (inputSize)
                    break;
                if (!run_)
                    return;
                usleep(1e3);
            }
            inputMutex_.lock();
            auto data = input_.front();
            inputMutex_.unlock();
            if (!(*chanIterator)->processSamples(data)) {
                cout << "PRN " << (*chanIterator)->prn() << " lost track" << endl;
                chanIterator = channels_.erase(chanIterator);
            } else
                ++chanIterator;
        }
        inputMutex_.lock();
        input_.pop();
        inputMutex_.unlock();
        if (userEstimateEcefTime_[3] != 0.0)
            userEstimateEcefTime_[3] += CA_CODE_TIME;
        if (++trackingPacketCount % 100 == 0) {
            solve();
        }
    }
}
