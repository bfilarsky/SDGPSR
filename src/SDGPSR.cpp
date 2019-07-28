#include "SDGPSR.h"
#include <iomanip>
#include <unistd.h>

const double SEARCH_WINDOW_BANDWIDTH = 10e3;
const double SEARCH_WINDOW_STEP_SIZE = 500.0;
const double WGS84_SEMI_MAJOR_AXIS   = 6378137.0;
const double WGS84_E_SQUARED         = 6.69437999014e-3;
const unsigned SEC_PER_WEEK          = 604800;
SDGPSR::SDGPSR(double fs, double clockOffset) : fs_(fs),
        fft_(fs_ * CA_CODE_TIME),
        run_(true),
        signalProcessor_(&SDGPSR::threadFunction, this) {
    clockOffset_ = clockOffset;
    userEstimateEcefTime_ = Vector4d(0.0, 0.0, 0.0, 0.0);
    navSolutionStarted_ = false;
    synced_ = true;

#ifdef DEBUG_FILES
    userEstimates_.open("userEstimates.bin", std::ofstream::binary);
#endif
}

SDGPSR::~SDGPSR() {
    run_ = false;
    signalProcessor_.join();
}

void SDGPSR::basebandSignal(const fftwVector &data) {
    std::lock_guard<std::mutex> lock(inputMutex_);
    input_.push(data);
}

void SDGPSR::basebandSignal(fftwVector &&data) {
    std::lock_guard<std::mutex> lock(inputMutex_);
    input_.push(std::move(data));
}

bool SDGPSR::synced(void) {
    std::lock_guard<std::mutex> lock(inputMutex_);
    return synced_ && !input_.size();
}

/*
 * Compute the circular correlation between searchData and basebandCode over a 1 ms interval corrCount times, and sum the magnitude
 * of each result. This allows the integration over periods that cross nav-bit boundaries at the expense of SNR
 */
std::vector<double> SDGPSR::nonCoherentCorrelator(std::vector<fftwVector> &searchData, fftwVector &basebandCode, unsigned corrCount) {
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

/*
 * Generate a replicated PRN code for the given prn and frequency offset
 */
void SDGPSR::basebandGenerator(unsigned prn, fftwVector &basebandCode, double freqOffset) {
    CaCode caCode(prn);
    basebandCode.resize(fs_ * CA_CODE_TIME);
    for (unsigned i = 0; i < basebandCode.size(); ++i) {
        double time = i / fs_;
        int chip = caCode[time * CHIP_RATE];
        double phase = 2.0 * M_PI * freqOffset * time;
        basebandCode[i] = std::complex<double>(chip * cos(phase), chip * sin(phase));
    }
}

/*
 * Conduct a 2-d search for the given prn using the frequency-domain searchData for all chip offsets and frequencies
 * between the start and stop at the given step size. The search will be conducted non-coherently using corrCount integrations
 */
SearchResult SDGPSR::search(std::vector<fftwVector> &searchData,
        unsigned prn,
        unsigned corrCount,
        double freqStart,
        double freqStop,
        double freqStep) {
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

#ifdef DEBUG_FILES
    std::string filename = "prn" + std::to_string(prn) + "searchWindow" + std::to_string((int) freqStep) + ".bin";
    std::ofstream output(filename, std::ofstream::binary);
    for (unsigned i = 0; i < searchWindow.size(); ++i)
        output.write((char*) &searchWindow[i][0], searchWindow[i].size() * sizeof(searchWindow[i][0]));
#endif

    //Calculate the statistics of the search pattern
    double mean = sum / (freqWindowSize * basebandCode.size());
    double accum = 0.0;
    for (unsigned freqWindow = 0; freqWindow < freqWindowSize; ++freqWindow)
        for (unsigned i = 0; i < searchWindow[freqWindow].size(); ++i)
            accum += (searchWindow[freqWindow][i] - mean) * (searchWindow[freqWindow][i] - mean);
    double stdDev = sqrt(accum / (freqWindowSize * basebandCode.size() - 1));
    double maxAboveMeanStdDevs = (max - mean) / stdDev;

    //If the largest value was at least SAT_FOUND_THRESH number of standard deviations above the mean value,
    //consider the satellite to be found at that point
    SearchResult searchResult;
    if (maxAboveMeanStdDevs > SAT_FOUND_THRESH) {
        searchResult.found = true;
        searchResult.baseBandFreq = peakFreq;
        searchResult.power = maxAboveMeanStdDevs;
        searchResult.sampleOffset = sampleOffset;
    }

    return searchResult;
}

/*
 * Attempt to solve for user time and position. The solution will only be calculated if at least 4
 * tracking channels have a full nav solution for the satellite (time, clock data, and ephemeris data)
 */
void SDGPSR::solve(void) {
    vector<pair<unsigned,Vector4d>> satPosAndTime;
    double maxTime = 0.0;
    for (auto &chan : channels_) {
        chan->sync();
        if (chan->state() == fullNav) {
            double time = chan->transmitTime();
            Vector3d position = chan->satellitePosition(time);
            //If we already have some idea of the nav solution, we can account for
            //the rotation of the earth between the time the satellite broadcast its signal
            //and we received it
            if (navSolutionStarted_) {
                double tof = userEstimateEcefTime_[3] - time;
                double earthRotation = -tof * OMEGA_EARTH;
                position[0] = position[0] * cos(earthRotation) - position[1] * sin(earthRotation);
                position[1] = position[0] * sin(earthRotation) + position[1] * cos(earthRotation);
            }
            Vector4d temp;
            temp[0] = position.x();
            temp[1] = position.y();
            temp[2] = position.z();
            temp[3] = time;
            if (time > maxTime)
                maxTime = time;
            satPosAndTime.emplace_back(chan->prn(),temp);
        }
    }

    if (satPosAndTime.size() >= 4) {
        if (!navSolutionStarted_){
            std::lock_guard<std::mutex> lock(userEstMutex_);
            userEstimateEcefTime_[3] = maxTime + .08;
            navSolutionStarted_ = true;
        }
        //Have over-determined system of equations residuals = H * userEstimateError, want to compute the least squares solution
        //to get userEstimate.
        //userEstimate += (H^T*H)^-1 * H^T*residuals.
        MatrixXd hMatrix(satPosAndTime.size(), 4);
        VectorXd residuals(satPosAndTime.size());
        for (unsigned i = 0; i < satPosAndTime.size(); ++i) {
            //Handle case where the week rollover has happened for the user, but the received signal is from the previous week
            double satelliteTime = satPosAndTime[i].second.w();
            if (navSolutionStarted_ && userEstimateEcefTime_[3] < SEC_PER_WEEK / 2 && satelliteTime > SEC_PER_WEEK / 2)
                satelliteTime -= SEC_PER_WEEK;
            //Compute range, pseudorange, and the residual (error between the two)
            double xRangeEst = satPosAndTime[i].second.x() - userEstimateEcefTime_.x();
            double yRangeEst = satPosAndTime[i].second.y() - userEstimateEcefTime_.y();
            double zRangeEst = satPosAndTime[i].second.z() - userEstimateEcefTime_.z();
            double rangeEstimate = sqrt(xRangeEst * xRangeEst + yRangeEst * yRangeEst + zRangeEst * zRangeEst);
            double pseudorangeEst = (userEstimateEcefTime_.w() - satelliteTime) * SPEED_OF_LIGHT_MPS;
            residuals(i) = rangeEstimate - pseudorangeEst;

#ifdef DEBUG_FILES
            if (!residualsOutput_.count(satPosAndTime[i].first))
                residualsOutput_.emplace(satPosAndTime[i].first, std::ofstream("residualsPrn" + std::to_string(satPosAndTime[i].first) + ".bin"));
            residualsOutput_[satPosAndTime[i].first].write((char*) &userEstimateEcefTime_[3], sizeof(userEstimateEcefTime_[3]));
            residualsOutput_[satPosAndTime[i].first].write((char*) &rangeEstimate, sizeof(rangeEstimate));
            residualsOutput_[satPosAndTime[i].first].write((char*) &pseudorangeEst, sizeof(pseudorangeEst));
#endif

            hMatrix(i, 0) = xRangeEst / rangeEstimate;
            hMatrix(i, 1) = yRangeEst / rangeEstimate;
            hMatrix(i, 2) = zRangeEst / rangeEstimate;
            hMatrix(i, 3) = 1.0;
        }
        //User cholesky decomposition solver to solve for the delta estimate
        Vector4d deltaEst = (hMatrix.transpose() * hMatrix).ldlt().solve(hMatrix.transpose() * residuals);
        std::lock_guard<std::mutex> lock(userEstMutex_);
        userEstimateEcefTime_[0] += deltaEst[0];
        userEstimateEcefTime_[1] += deltaEst[1];
        userEstimateEcefTime_[2] += deltaEst[2];
        userEstimateEcefTime_[3] += deltaEst[3] / SPEED_OF_LIGHT_MPS;

#ifdef DEBUG_FILES
        userEstimates_.write((char*) &userEstimateEcefTime_, sizeof(userEstimateEcefTime_));
#endif
    }
}

void SDGPSR::threadFunction() {
    //FFT the first CORR_COUNT packets in order to conduct search
    const unsigned CORR_COUNT = 128;
    std::vector<fftwVector> searchData(CORR_COUNT);
    for (unsigned i = 0; i < CORR_COUNT; ++i) {
        //Wait for available data
        while (1) {
            inputMutex_.lock();
            size_t inputSize = input_.size();
            inputMutex_.unlock();
            if (inputSize)
                break;
            if (!run_)
                return;
            synced_ = true;
            usleep(1e3);
        }
        synced_ = false;
        //Once data is available, fft it and output it the search buffer
        inputMutex_.lock();
        searchData[i].resize(input_.front().size());
        fft_.forward(&input_.front()[0], &searchData[i][0]);
        input_.pop();
        inputMutex_.unlock();
    }

    //Using the FFT'd data, search for all PRNs
    for (unsigned prn = 1; prn <= 32; ++prn) {
        SearchResult searchResult = search(searchData,
                prn,
                CORR_COUNT,
                clockOffset_ - SEARCH_WINDOW_BANDWIDTH / 2.0,
                clockOffset_ + SEARCH_WINDOW_BANDWIDTH / 2.0,
                SEARCH_WINDOW_STEP_SIZE);
        if (searchResult.found) {
            std::lock_guard<std::mutex> lock(channelMutex_);
            channels_.emplace_back(new SignalTracker(fs_, prn, searchResult));
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
                synced_ = true;
                usleep(1e3);
            }
            synced_ = false;
            inputMutex_.lock();
            State state = (*chanIterator)->processSamples(input_.front());
            inputMutex_.unlock();
            if (state == lossOfLock) {
                std::lock_guard<std::mutex> lock(channelMutex_);
                chanIterator = channels_.erase(chanIterator);
            } else
                ++chanIterator;
        }
        inputMutex_.lock();
        input_.pop();
        inputMutex_.unlock();
        if (navSolutionStarted_){
            std::lock_guard<std::mutex> lock(userEstMutex_);
            userEstimateEcefTime_[3] = fmod(userEstimateEcefTime_[3] + CA_CODE_TIME, SEC_PER_WEEK);
        }
        if (++trackingPacketCount % 100 == 0) {
            solve();
        }
    }
}

Vector3d SDGPSR::positionECEF(void){
    std::lock_guard<std::mutex> lock(userEstMutex_);
    return Vector3d(userEstimateEcefTime_.x(), userEstimateEcefTime_.y(), userEstimateEcefTime_.z());
}

Vector3d SDGPSR::positionLLA(void){
    Vector3d posECEF          = positionECEF();
    if (posECEF == Vector3d(0.0,0.0,0.0))
        return posECEF;
    double longitude          = atan2(posECEF.y(), posECEF.x());
    double xyPlaneRadius      = sqrt(posECEF.x() * posECEF.x() + posECEF.y() * posECEF.y());
    double geocentricLatitude = atan2(xyPlaneRadius, posECEF.z());

    double geodeticLatitude   = geocentricLatitude;
    double height             = 0.0;

    while (true){
        double oldGeodeticLatitude = geodeticLatitude;
        double sinGeodetic         = sin(geodeticLatitude);
        double radiusOfCurvature   = WGS84_SEMI_MAJOR_AXIS / sqrt(1.0 - WGS84_E_SQUARED * sinGeodetic * sinGeodetic);
        height                     = xyPlaneRadius/cos(geodeticLatitude) - radiusOfCurvature;
        geodeticLatitude           = atan(posECEF.z() / (xyPlaneRadius * (1 - WGS84_E_SQUARED * radiusOfCurvature / (radiusOfCurvature + height))));
        //Want values to converge to 6 decimal places, or < 6 inches
        if (fabs(oldGeodeticLatitude - geodeticLatitude) < 1e-6)
            break;
    }

    return Vector3d(geodeticLatitude * 180.0 / M_PI, longitude * 180.0 / M_PI, height);
}

double SDGPSR::timeOfWeek(void){
    std::lock_guard<std::mutex> lock(userEstMutex_);
    return userEstimateEcefTime_.w();
}

std::vector<std::pair<unsigned, State>> SDGPSR::trackingStatus(){
    std::vector<std::pair<unsigned, State>> status;
    std::lock_guard<std::mutex> lock(channelMutex_);
    for (auto &channel : channels_)
        status.emplace_back(channel->prn(), channel->state());
    return status;
}

bool SDGPSR::navSolution(void){
    return navSolutionStarted_;
}
