#include <iostream>
#include <memory>
#include <chrono>

#include "AudioDevice.hpp"
#include "SpectrumAnalyzer.hpp"

constexpr int SAMPLE_RATE = 48000;
constexpr int AUDIO_BLOCK_SIZE = 1024;
constexpr int SPECTRUM_BLOCK_SIZE = 4096;

int main() {
    auto audioDevice = std::make_shared<AudioDevice>(AudioDevice::DEFAULT_DEVICE, SAMPLE_RATE, AUDIO_BLOCK_SIZE);
    auto spectrumAnalyzer = std::make_shared<SpectrumAnalyzer>(audioDevice, SPECTRUM_BLOCK_SIZE);

    audioDevice->startStream();

    while(true) {
        std::cout << "Test" << std::endl;

        std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }
}