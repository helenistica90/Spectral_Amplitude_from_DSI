# Spectral_Amplitude_from_DSI
MATLAB script for computing DSI (DisSimilarity Index) from Rayleigh wave Ellipticity. Includes time-series resampling, periodicity detection via FFT, and extraction of 12-hour &amp; 24-hour spectral components. Data required: timestamps, ellipticity matrix, and frequency vector.

# DSI Analysis from Rayleigh Wave Ellipticity

## 📌 Overview
This repository contains the MATLAB script for computing the DisSimilarity Index (DSI) from Rayleigh wave ellipticity. The methodology is designed to get the time-lapse spectral amplitude, derived from DSI series as proxy of velocity variations, at 1 cycle/day and 2 cycles/day periodicities.
## 📊 Features
- Computes DSI time series from Rayleigh wave ellipticity data
- Performs time-series resampling and smoothing
- Applies Fourier analysis to detect periodic signals
- Extracts 12-hour and 24-hour periodic components

## 📂 Required Data
- **dateX**: Time vector representing the timestamps
- **MatFreq**: Matrix storing ellipticity values per frequency band
- **Frecusmedia**: Frequency vector associated with MatFreq

## 🛠 Usage
1. Load the provided `.mat` dataset.
2. Run the `DSI_and_SpectralAmplitude.m` script.
3. Results are saved in `DSI_results.mat` and `Ampli_12_results.mat`.

## 🔬 Reference
If you use this code for research, please cite:  
Seivane H. (2025). DSI Analysis from Rayleigh Wave Ellipticity. Available at: https://github.com/helenistica90/Spectral_Amplitude_from_DSI

## 📧 Contact
For questions or collaborations, contact helenaseiv@outlook.com


