# ECG Signal Integrity: Design and Verification of a Digital Biomedical Signal Processing Pipeline

**Author:** Pirrijean  
**Framework:** [ECG-Integrity-Lab](https://github.com/Pirrijean/ECG-Integrity-Lab) (`v0.1.0`)  

---

## Abstract / Executive Summary

How reliable is a classical electrocardiogram (ECG) QRS detection algorithm when the physiological input is degraded by power-line interference, respiratory baseline drift, and broadband instrumentation noise? To investigate this question rigorously without requiring clinical acquisition hardware, **ECG-Integrity-Lab** was developed in Python as an open biomedical signal conditioning and verification framework using expert-annotated recordings from the **MIT-BIH Arrhythmia Database** ($f_s = 360\text{ Hz}$, 11-bit ADC resolution). Across a five-record full-duration pilot encompassing **10,789 expert-reviewed cardiac beats** (64,734 beat evaluations across six controlled noise scenarios), a zero-phase 4th-order Butterworth bandpass ($5\text{–}20\text{ Hz}$) paired with a Pan–Tompkins-inspired adaptive detector completely suppressed out-of-band $0.25\text{ Hz}$ baseline wander ($0.20\text{ mV}$) and $60\text{ Hz}$ mains interference ($0.05\text{ mV}$) with zero added false detections. However, broadband Gaussian noise ($20\text{ dB} \rightarrow 5\text{ dB}$ SNR) and clinical in-band electrode motion artifacts (`Record 105`) or multiform premature ventricular contractions (`Record 106`) revealed critical degradation mechanisms: temporal R-peak localization error ($\text{MAE}_t$) rose by up to $+190\%$ (from $1.23\text{ ms}$ to $3.56\text{ ms}$ on `Record 105`, and up to $4.27\text{ ms}$ on `Record 106`, surpassing the $2.78\text{ ms}$ sampling interval), while missed ectopic beats tripled at $5\text{ dB}$ SNR.

---

## 1. Problem Statement: Signal Integrity in Automated ECG Analysis

Biomedical surface potentials acquired at the skin–electrode interface have millivolt-scale amplitudes ($0.5\text{–}3\text{ mV}$) that are highly vulnerable to coupled electrical and physiological disturbances:

- **Baseline wander ($< 0.5\text{ Hz}$):** Induced by thoracic respiration and skin–electrode half-cell potential shifts.
- **Power-line interference ($50/60\text{ Hz}$):** Capacitive and inductive mains coupling alongside finite amplifier CMRR.
- **Broadband electromyographic (EMG) & thermal noise:** Spectral overlap directly inside the $5\text{–}40\text{ Hz}$ QRS energy band.
- **Electrode motion artifacts ($1\text{–}15\text{ Hz}$):** Transient mechanical impedance changes mimicking QRS complexes in amplitude and slope.

Evaluating an algorithm exclusively on clean sinus rhythm risks overstating performance. Engineering verification requires controlled perturbation testing paired with quantitative beat-matching and temporal localization metrics.

---

## 2. Methodology: Dataset Audit & Calibration QA

Recordings were imported from PhysioNet's **MIT-BIH Arrhythmia Database** (`mitdb`, ODC-By v1.0) via WFDB Python. Prior to filtering, an automated QA gate (`dataset.py`) audited each header:

- Confirmed $f_s = 360.0\text{ Hz}$ ($\Delta t = 2.7778\text{ ms}$), physical calibration units (`mV`), `200.0 counts/mV` ADC gain, and `1024` digital baseline counts—never treating uncalibrated ADC integers as millivolts.
- Verified absence of `NaN`/`Inf` samples and enforced physiological amplitude bounds ($|x| \le 20\text{ mV}$).
- Filtered `.atr` annotations strictly to standard PhysioNet QRS beat symbols, excluding non-beat rhythm/quality markers (e.g., `119` artifact/noise markers in `Record 105` and `71` rhythm transitions in `Record 106`).

---

## 3. Pipeline Architecture: Conditioning, Perturbation & QRS Detection

1. **Controlled Perturbation Generator (`signal_processing.py`):** Applies deterministic $0.25\text{ Hz}$ wander ($0.20\text{ mV}$), $60\text{ Hz}$ mains ($0.05\text{ mV}$), and seeded zero-mean Gaussian noise (`seed=42`) scaled to AC signal variance at $20$, $10$, and $5\text{ dB}$ SNR.
2. **QRS-Oriented Conditioning:** 4th-order Butterworth bandpass ($5\text{–}20\text{ Hz}$) in second-order sections (`sosfiltfilt`), achieving zero phase distortion while attenuating out-of-band drift and mains hum. *(Explicitly scoped for QRS localization only; unsuitable for diagnostic ST/P/T morphology).*
3. **Pan–Tompkins-Inspired Detector (`detector.py`):** Centered differentiation $d[n]$, squaring $e[n] = d[n]^2$, $120\text{ ms}$ moving-window integration $M[n]$, dual adaptive signal/noise thresholds (`SPKI`/`NPKI`) with $1.66 \times \overline{RR}$ searchback, $250\text{ ms}$ refractory period, and $\pm 90\text{ ms}$ R-peak relocation.

---

## 4. Software Metrological Verification (`digital_phantom.py`)

Before analyzing human ECGs, amplitude and time-base routines were verified against a synthetic $1.0\text{ mV}$, $60.0\text{ bpm}$, $100\text{ ms}$ rectangular pulse train ($10\text{ s}$ at $360\text{ Hz}$):

- Nominal pulse count (`10`), amplitude error ($0.0\%$), and rate error ($0.0\%$) were confirmed exactly.
- Injecting an artificial $+5\%$ scale gain ($1.05\text{ mV}$) yielded an exact $+5.0\%$ measured amplitude error ($E_A$), verifying software sensitivity to calibration drift (while clearly distinguishing digital software verification from SI-traceable physical hardware calibration).

---

## 5. Quantitative Evaluation Results

Using strict **one-to-one beat matching** within a $\pm 150\text{ ms}$ tolerance window (`evaluation.py`), the five-record 30-minute pilot (`100`, `101`, `103`, `105`, `106`) demonstrated:

- **Complete rejection of out-of-band sinusoidal artifacts:** Both $0.25\text{ Hz}$ baseline wander and $60\text{ Hz}$ mains interference produced **zero additional False Positives** across all 10,789 beats.
- **Temporal jitter precedes beat-loss under broadband noise:** On clean `Record 100` ($n = 2,265$ beats), F1-score remained $100.00\%$ even at $5\text{ dB}$ SNR, yet $\text{MAE}_t$ rose monotonically from $1.50\text{ ms}$ to $1.96\text{ ms}$ ($+30\%$).
- **Vulnerability to in-band motion artifacts & wide PVCs:**
  - `Record 105` (high-grade motion noise) exhibited `33 FP` and `8 FN` at baseline ($\text{PPV} = 98.73\%$, $F_1 = 99.20\%$), with $\text{MAE}_t$ surging from $1.23\text{ ms}$ to **$3.56\text{ ms}$** at $5\text{ dB}$ SNR.
  - `Record 106` ($518$ multiform PVCs) suffered a $3\times$ increase in missed beats (`2 FN` $\rightarrow$ `6 FN`) and `4 FP` at $5\text{ dB}$ SNR, with $\text{MAE}_t$ reaching **$4.27\text{ ms}$** ($> 1.5 \times \Delta t$) due to reduced high-frequency $dV/dt$ slope in wide ventricular ectopic complexes.

---

## 6. Conclusions & Portfolio Summary Post (Ready for Technical Sharing)

```text
ECG Signal Integrity: Building a Digital Verification Framework for Biomedical Signals

How reliable is an ECG processing algorithm when the input signal is affected by electrical interference, physiological noise, and arrhythmia morphology changes?

I developed ECG-Integrity-Lab, an open Python framework to investigate this question quantitatively using expert-annotated recordings from PhysioNet's MIT-BIH Arrhythmia Database (360 Hz, 11-bit ADC).

Across a 5-record full-duration pilot (10,789 expert-reviewed beats; 64,734 beat evaluations across 6 controlled perturbation scenarios), the framework integrates:

• Biomedical instrumentation & QA: Automated header verification of sampling rate, physical units (mV), ADC gain (200 counts/mV), digital baseline, and exclusion of non-beat rhythm/noise markers.
• Software metrology: A deterministic 1.0 mV / 60 bpm digital calibration phantom (0.0% nominal scale/rate error; verified +5.0% gain-error detection) plus automated unit tests.
• Digital signal conditioning: Zero-phase 4th-order Butterworth SOS bandpass (5–20 Hz) paired with controlled baseline wander (0.25 Hz), 60 Hz mains hum, and seeded AWGN (20, 10, and 5 dB SNR).
• Quantitative failure analysis: Strict 1-to-1 beat matching (±150 ms) revealed zero false positives from out-of-band 0.25 Hz / 60 Hz interference, whereas broadband 5 dB noise and clinical motion/ectopic morphologies (Records 105 & 106) increased R-peak temporal localization error (MAE_t) up to 3.56–4.27 ms (>1 sample period) and tripled missed multiform PVCs.

This is an exploratory biomedical engineering research and software-verification project, not a clinically validated diagnostic system or medical device.

Repository, protocol & reproducible results: https://github.com/Pirrijean/ECG-Integrity-Lab

#BiomedicalEngineering #BiomedicalInstrumentation #ECG #SignalProcessing #MedicalDevices #QualityEngineering #Python #Metrology #RAndD
```
