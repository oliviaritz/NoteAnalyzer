#include "SpectrumAnalyzer.hpp"

#include <iostream>

using namespace std;

Spectrum::Spectrum(unsigned int binCount, unsigned int sampleRate)
    : m_bins(binCount)
    , m_sampleRate{ sampleRate } {
}

float Spectrum::getMinFrequency() const {
    return 0;
}

float Spectrum::getMaxFrequency() const {
    return m_sampleRate / 2;
}

float Spectrum::getFrequencyOfBin(unsigned int bin) const {
    return bin * m_sampleRate / (m_bins.size() * 2.f);
}

unsigned int Spectrum::getBinOfFrequency(float freq) const {
    return freq * (m_bins.size() * 2.f) / m_sampleRate;
}

float Spectrum::getPeakBinFrequency(float minFreq, float maxFreq) const {
    auto startBin = getBinOfFrequency(minFreq);
    auto endBin = getBinOfFrequency(maxFreq);

    auto maxIndex = std::max_element(m_bins.begin() + startBin, m_bins.begin() + endBin) - m_bins.begin();

    return getFrequencyOfBin(maxIndex);
}

SpectrumAnalyzer::SpectrumAnalyzer(std::shared_ptr<AudioDevice>& _audioDevice,
	unsigned int _blockSize)
	:	workUnit(std::make_unique<boost::asio::io_service::work>(ioService))
	,	audioDevice(_audioDevice)
    , leftSpectrum{ _blockSize / 2, audioDevice->getSampleRate() }
    , rightSpectrum{ _blockSize / 2, audioDevice->getSampleRate() }
	,	chunkSize{audioDevice->getBlockSize()}
	,	blockSize{_blockSize} {

	//Initialize FFTW stuff
	//fftIn, fftOut, fftPlan
	fftIn = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * blockSize);
	fftOut = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * blockSize);

	//Compute FFT plans
	fftPlan = fftw_plan_dft_1d(blockSize, fftIn, fftOut,
		FFTW_FORWARD, FFTW_MEASURE);

	//Generate FFT window function
	generateWindow();

	//Initialize audio buffers
	leftBuffer.resize(blockSize);
	rightBuffer.resize(blockSize);

	//Launch thread
	asyncThread = std::thread(std::bind(&SpectrumAnalyzer::threadRoutine,
		this));

	//Register audio callback
	callbackID = audioDevice->addCallback(std::bind(&SpectrumAnalyzer::cbAudio, this,
		std::placeholders::_1, std::placeholders::_2));
}

SpectrumAnalyzer::~SpectrumAnalyzer() {
	//Remove audio callback
	audioDevice->removeCallback(callbackID);

	//Shutdown threads
	workUnit.reset();

	asyncThread.join();

	//Cleanup fftw stuff
	fftw_destroy_plan(fftPlan);
	fftw_free(fftIn);
	fftw_free(fftOut);
}

std::shared_ptr<AudioDevice> SpectrumAnalyzer::getAudioDevice() {
	return audioDevice;
}

Spectrum SpectrumAnalyzer::getLeftSpectrum() const {
	std::unique_lock<std::mutex> spectrumLock(spectrumMutex);

	return leftSpectrum;
}

Spectrum SpectrumAnalyzer::getRightSpectrum() const {
	std::unique_lock<std::mutex> spectrumLock(spectrumMutex);

	return rightSpectrum;
}

void SpectrumAnalyzer::cbAudio(const int16_t* left, const int16_t* right) {
	std::unique_lock<std::mutex> bufferLock(bufferMutex);

	//Shift the samples forward by 1 chunk size
	std::memmove(leftBuffer.data(), &leftBuffer[chunkSize],
		sizeof(int16_t) * (blockSize - chunkSize));
	std::memmove(rightBuffer.data(), &rightBuffer[chunkSize],
		sizeof(int16_t) * (blockSize - chunkSize));

	//Copy the new chunk samples into the end of the block
	std::memcpy(&leftBuffer[blockSize - chunkSize], left,
		sizeof(int16_t) * chunkSize);
	std::memcpy(&rightBuffer[blockSize - chunkSize], right,
		sizeof(int16_t) * chunkSize);

	//Post the fft routine to the async thread
	ioService.post(std::bind(&SpectrumAnalyzer::fftRoutine, this,
		leftBuffer, rightBuffer));
}

void SpectrumAnalyzer::threadRoutine() {
	//Run work from ioService
	ioService.run();
}

void SpectrumAnalyzer::fftRoutine(std::vector<int16_t> left,
	std::vector<int16_t> right) {

	double sampleRate = audioDevice->getSampleRate();
	
	{
		std::unique_lock<std::mutex> bufferLock(bufferMutex);
		
		//Fill FFT input buffer
		for(unsigned int i = 0; i < blockSize; i++) {
			fftIn[i][0] = fftWindow[i] *
				((double)left[i] / INT16_MAX / blockSize);
			fftIn[i][1] = 0.; //Imaginary
		}
	}

	//Do FFT on left samples
	fftw_execute(fftPlan);

	//Fill left spectrum with new FFT data
	for(unsigned int i = 0; i < blockSize/2; ++i) {
        leftSpectrum.m_bins[i] = std::sqrt(sqr(fftOut[i][0]) +
            sqr(fftOut[i][1]));
	}
	
	//Now do right FFT

	{
		std::unique_lock<std::mutex> bufferLock(bufferMutex);

		//Copy real audio data into complex fft input array and scale to [-1., 1.]
		for(unsigned int i = 0; i < blockSize; i++) {
			fftIn[i][0] = fftWindow[i] * ((double)right[i] / INT16_MAX / blockSize); //Real
			fftIn[i][1] = 0.; //Imaginary
		}
	}

	//Do FFT on right samples
	fftw_execute(fftPlan);

	//Fill right spectrum with new FFT data
    for (unsigned int i = 0; i < blockSize / 2; ++i) {
        leftSpectrum.m_bins[i] = std::sqrt(sqr(fftOut[i][0]) +
            sqr(fftOut[i][1]));
    }
}

void SpectrumAnalyzer::generateWindow() {
	//Hanning window
	fftWindow.resize(blockSize);

	for(unsigned int i = 0; i < blockSize; i++) {
		fftWindow[i] = 0.5 * (1. - std::cos((2*3.141592654*i)/(blockSize - 1)));
	}
}

//Helper function
double SpectrumAnalyzer::sqr(const double x) {
	return x*x;
}
