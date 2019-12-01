#pragma once

#include <memory>
#include <thread>
#include <functional>
#include <cstdint>
#include <vector>

#include <boost/asio.hpp>

#include <fftw3.h>

#include "AudioDevice.hpp"


class Spectrum
{
public:
	Spectrum(unsigned int binCount, unsigned int sampleRate);
	
	float getMinFrequency() const;
	float getMaxFrequency() const;

    float getFrequencyOfBin(unsigned int bin) const;
    unsigned int getBinOfFrequency(float freq) const;

	float getPeakBinFrequency(float minFreq, float maxFreq) const;

private:
    friend class SpectrumAnalyzer;

	std::vector<float> m_bins;
	unsigned int m_sampleRate;
};


class SpectrumAnalyzer
{
public:
	SpectrumAnalyzer(std::shared_ptr<AudioDevice>& audioDevice, unsigned int blockSize);
	~SpectrumAnalyzer();


	std::shared_ptr<AudioDevice> getAudioDevice();

	Spectrum getLeftSpectrum() const;
	Spectrum getRightSpectrum() const;

private:
	void threadRoutine();
	void cbAudio(const int16_t* left, const int16_t* right);
	void fftRoutine(std::vector<int16_t>, std::vector<int16_t>);
	void generateWindow();

	static double sqr(const double x);

	//Thread stuff
	boost::asio::io_service ioService;
	std::unique_ptr<boost::asio::io_service::work> workUnit;
	std::thread asyncThread;

    std::shared_ptr<AudioDevice> audioDevice;

	Spectrum leftSpectrum, rightSpectrum;
	mutable std::mutex spectrumMutex;

	//Audio sample buffers
	std::vector<int16_t> leftBuffer, rightBuffer;
	std::mutex bufferMutex;

	//FFT stuff
	fftw_complex *fftIn, *fftOut;
	fftw_plan fftPlan;
	std::vector<double> fftWindow;

	//Audio device stuff
	unsigned int callbackID;
	unsigned int chunkSize; //Size of buffer from audio device
	unsigned int blockSize;	//Size of buffer sent through fft
};
