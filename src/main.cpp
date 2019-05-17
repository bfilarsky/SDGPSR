#include <iostream>
#include <fstream>
#include <unistd.h>
#include "libhackrf/hackrf.h"
#include "SDGPSR.h"
#include "LowPassFilter.h"

using std::cout;
using std::endl;

uint8_t *ptr = NULL;
unsigned position = 0;

const unsigned SAMPLE_RATE = 8e6;
const unsigned DECIMATION = 1;
const unsigned PACKET_SIZE = SAMPLE_RATE * CA_CODE_TIME;
const unsigned TIME = 200;
const unsigned BUFFER_SIZE = TIME * SAMPLE_RATE;
const unsigned THROWAWAYS = 5*SAMPLE_RATE;
unsigned throwawayCount = 0;

SDGPSR sdgpsr(SAMPLE_RATE / DECIMATION);


int rx_callback(hackrf_transfer* transfer){
	if (throwawayCount < THROWAWAYS){
		throwawayCount += transfer->buffer_length;
		return 0;
	}
	size_t copyLength = std::min((int)(BUFFER_SIZE - position), transfer->buffer_length);
	memcpy(&ptr[position], transfer->buffer, copyLength);
	position += copyLength;
	return 0;
	//return sdgpsr.basebandSignal(transfer->buffer, transfer->buffer_length);
}

void runHackRf(void){
	hackrf_init();
	hackrf_device_list_t *dlPtr = hackrf_device_list();
	hackrf_device* device = NULL;
	if (hackrf_device_list_open(dlPtr, 0, &device) != HACKRF_SUCCESS){
		cout << "Could not open device" << endl;
		exit(1);
	}
	hackrf_device_list_free(dlPtr);
	hackrf_set_freq(device, GPS_L1_HZ);
	hackrf_set_sample_rate(device, SAMPLE_RATE);
	uint32_t bandwidth_hz = 0;
	hackrf_compute_baseband_filter_bw(bandwidth_hz);
	hackrf_set_baseband_filter_bandwidth(device, bandwidth_hz);
	hackrf_set_vga_gain(device, 46);  //0-62, 2 dB steps
	hackrf_set_lna_gain(device, 40);  //0-40, 8 dB steps
	hackrf_set_antenna_enable(device, 1);
	hackrf_set_amp_enable(device, 0);
	ptr = (uint8_t *)malloc(BUFFER_SIZE);
	if (!ptr)
		return;
	hackrf_start_rx(device, rx_callback, NULL);
	while(position < BUFFER_SIZE);
	hackrf_stop_rx(device);
	hackrf_close(device);
	std::ofstream output("/home/bfilarsky/eclipse-workspace/SDGPSR/Debug/output.bin", std::ofstream::binary);
	output.write((char*)ptr, BUFFER_SIZE);
	output.close();
	free(ptr);
}

void runFromFile(void){
	std::ifstream input("/home/bfilarsky/eclipse-workspace/SDGPSR/Debug/output.bin", std::ifstream::binary | std::ios::ate);
	//std::ifstream input("/home/bfilarsky/output.bin", std::ifstream::binary | std::ios::ate);
	if (input.is_open()){
		size_t sizeBytes = input.tellg();
		size_t size = sizeBytes / sizeof(std::complex<int8_t>);
		input.seekg(0);
		//double lpf = 0.1;
		//LowPassFilter real(lpf);
		//LowPassFilter imag(lpf);
		unsigned packetsToRead = size / PACKET_SIZE;
		for (unsigned i = 0; i < packetsToRead && !input.eof(); ++i){
			std::vector<std::complex<int8_t>> intData(PACKET_SIZE);
			input.read((char*)&intData[0], PACKET_SIZE * sizeof(intData[0]));
			fftwVector data(PACKET_SIZE / DECIMATION);
			for (unsigned j = 0; j < PACKET_SIZE / DECIMATION; ++j)
				data[j] = std::complex<double>(intData[j * DECIMATION].real(), intData[j * DECIMATION].imag());
			//cout << (int)intData[0].real() << ',' << (int)intData[0].imag() << endl;
			sdgpsr.basebandSignal(data);
		}
		input.close();

		sdgpsr.signalProcessing();
	}
	else {
		cout << "Could not open input file" << endl;
	}
}

int main(){
	//runHackRf();
	runFromFile();
	return 0;
}
