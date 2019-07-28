#include <iostream>
#include <fstream>
#include <unistd.h>
#include "src/SDGPSR.h"

using std::cout;
using std::endl;

const unsigned SAMPLE_RATE = 2e6;
const unsigned PACKET_SIZE = SAMPLE_RATE * CA_CODE_TIME;

std::vector<std::string> stateMap = {{"Loss of Lock"}, {"Closing Carrier Frequency-Locked Loop"}, {"Closing Carrier Phase-Locked Loop"}, {"Finding Nav-Bit Edge"}, {"Full track - Downloading required Nav Data"}, {"Full Nav Solution"}, };

int main(){
	//The SDR used has a clock error of approximately -10 KHz
	SDGPSR sdgpsr(SAMPLE_RATE, -10e3);
	std::ifstream input("sampleData/iqData.dat", std::ifstream::binary | std::ios::ate);
	if (input.is_open()){
		cout << "Loading data from file" << endl;
		size_t sizeBytes = input.tellg();
		size_t size = sizeBytes / sizeof(std::complex<int8_t>);
		input.seekg(0);
		unsigned packetsToRead = size / PACKET_SIZE;
		for (unsigned i = 0; i < packetsToRead && !input.eof(); ++i){
			std::vector<std::complex<int8_t>> intData(PACKET_SIZE);
			input.read((char*)&intData[0], PACKET_SIZE * sizeof(intData[0]));
			fftwVector data(PACKET_SIZE);
			for (unsigned j = 0; j < PACKET_SIZE; ++j)
				data[j] = std::complex<double>(intData[j].real(), intData[j].imag());
			sdgpsr.basebandSignal(std::move(data));
		}
		input.close();

		while(!sdgpsr.synced()){
			if (!sdgpsr.navSolution()){
				cout << "Acquiring Nav Solution. Satellite tracking states: " << endl;
				auto status = sdgpsr.trackingStatus();
				for (auto &prnStatus : status)
					cout << "PRN " << prnStatus.first << " is in state " << prnStatus.second << '(' << stateMap[prnStatus.second] << ')' << endl;
				cout << endl;
			}
			else {
				cout << "User Time of Week: " << sdgpsr.timeOfWeek() << endl;
		    	cout << "User Latitude (deg): "<< std::setprecision(10) << sdgpsr.positionLLA().x() << endl;
				cout << "User Longitude (deg): "<< std::setprecision(10) << sdgpsr.positionLLA().y() << endl;
				cout << "User Altitude (m): " << std::setprecision(10) << sdgpsr.positionLLA().z() << endl;
				cout << endl;
			}
		    usleep(1e6);
		}

	}
	else {
		cout << "Could not open input file" << endl;
	}
	return 0;
}
