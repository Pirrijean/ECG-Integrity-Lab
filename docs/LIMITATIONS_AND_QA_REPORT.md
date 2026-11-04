# Quantitative QA Report, Failure Mode Analysis & Engineering Limitations

**Project:** ECG-Integrity-Lab (`v0.1.0`)  
**Scope:** Empirical Evaluation of MIT-BIH Record `100` (5-min initial test) and 5-Record Full-Duration Pilot (`100`, `101`, `103`, `105`, `106` — 1800 s / 30 min each, Channel 0 `MLII`)  

---

## 1. Software & Metrological Pre-Flight Verification Results

Prior to processing physiological recordings, the software arithmetic and scaling pipeline were verified using `digital_phantom.py` (`python examples/run_digital_phantom.py`) and the automated `pytest` suite (`7/7` unit & negative fault-injection tests passing):

| Verification Test | Expected Value | Observed Value | Relative Error | Status |
| :--- | :---: | :---: | :---: | :---: |
| Nominal Pulse Count ($10\text{ s}$ @ $60\text{ bpm}$) | `10` | `10` | `0.00%` | PASS |
| Nominal Amplitude ($1.0\text{ mV}$) | `1.0000 mV` | `1.0000 mV` | `0.00%` | PASS |
| Nominal Heart Rate ($60.0\text{ bpm}$, $f_s = 360\text{ Hz}$) | `60.0000 bpm` | `60.0000 bpm` | `0.00%` | PASS |
| Scale Sensitivity Check ($+5\%$ Gain, $1.05\text{ mV}$) | `+5.00%` | `+5.00%` | `0.00%` | PASS |
| Fault Injection (`fs <= 0`, `NaN`, invalid channel, raw ADC counts $>20\text{ mV}$) | `ValueError` | `ValueError` | — | PASS |

---

## 2. Header & Data Integrity Audit Across Pilot Records (Step 9)

All five recordings satisfied physical calibration checks (`fs = 360.0 Hz`, $\Delta t = 2.7778\text{ ms}$, physical units `mV`, ADC gain `200.0 counts/mV`, 11-bit resolution, digital baseline `1024 counts`, `0` non-finite samples). Filtering `.atr` files strictly to true QRS beat symbols prevented rhythm/noise markers (`+`, `~`, `|`) from corrupting ground-truth counts:

| Record | Duration (s) | Samples | True Reference Beats (`n_ref`) | Excluded Non-Beat Markers | Beat Symbol Breakdown | Clinical / Technical Profile (MIT-BIH Header) |
| :---: | :---: | :---: | :---: | :---: | :--- | :--- |
| **`100`** (5 min) | 300.0 | 108,000 | 371 | 1 | `N: 367, A: 4` | Clean normal sinus rhythm, rare APCs |
| **`100`** (30 min) | 1800.0 | 648,000 | 2,265 | 1 | `N: 2231, A: 33, V: 1` | High-SNR sinus reference |
| **`101`** (30 min) | 1800.0 | 648,000 | 1,860 | 11 | `N: 1855, A: 3, Q: 2` | Sinus rhythm with baseline drift & isolated artifact bursts |
| **`103`** (30 min) | 1800.0 | 648,000 | 2,078 | 12 | `N: 2076, A: 2` | Clean sinus rhythm with subtle amplitude transitions |
| **`105`** (30 min) | 1800.0 | 648,000 | 2,565 | **119** | `N: 2519, V: 41, Q: 5` | **Predominant high-grade electrode motion noise & uniform PVCs** |
| **`106`** (30 min) | 1800.0 | 648,000 | 2,021 | **71** | `N: 1503, V: 518` | **Frequent multiform PVCs & ventricular bigeminy/trigeminy** |

> [!IMPORTANT]
> Had the 119 non-beat signal-quality/artifact markers in Record `105` or the 71 rhythm-transition markers in Record `106` been treated naively as QRS complexes, Sensitivity would have been artificially penalized by up to $4.4\%$ solely due to annotation misinterpretation rather than detector error.

---

## 3. Empirical Results: Initial 5-Minute Verification (Record `100`, `MLII`)

| Scenario | Condition | TP | FP | FN | Sensitivity (%) | PPV (%) | F1-Score (%) | $\text{MAE}_t$ (ms) |
| :--- | :--- | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
| **Baseline** | Original calibrated signal | 371 | 0 | 0 | 100.00 | 100.00 | 100.00 | 1.17 |
| **Baseline wander** | $0.20\text{ mV}$, $0.25\text{ Hz}$ | 371 | 0 | 0 | 100.00 | 100.00 | 100.00 | 1.23 |
| **Mains interference** | $0.05\text{ mV}$, $60.0\text{ Hz}$ | 371 | 0 | 0 | 100.00 | 100.00 | 100.00 | 1.26 |
| **White noise (20 dB)** | AWGN $20\text{ dB}$ SNR (`seed=42`) | 371 | 0 | 0 | 100.00 | 100.00 | 100.00 | 1.23 |
| **White noise (10 dB)** | AWGN $10\text{ dB}$ SNR (`seed=42`) | 371 | 0 | 0 | 100.00 | 100.00 | 100.00 | 1.50 |
| **White noise (5 dB)** | AWGN $5\text{ dB}$ SNR (`seed=42`) | 371 | 0 | 0 | 100.00 | 100.00 | 100.00 | 1.70 |

---

## 4. Empirical Results: 5-Record Full-Duration Pilot (10,789 Reference Beats)

Summary of F1-score, error counts (`FP` / `FN`), and temporal localization precision ($\text{MAE}_t$) across `results/pilot/pilot_metrics.csv`:

| Record | Baseline (`FP`/`FN`, F1, $\text{MAE}_t$) | Wander $0.25\text{ Hz}$ (`FP`/`FN`, F1, $\text{MAE}_t$) | Mains $60\text{ Hz}$ (`FP`/`FN`, F1, $\text{MAE}_t$) | AWGN $20\text{ dB}$ (`FP`/`FN`, F1, $\text{MAE}_t$) | AWGN $10\text{ dB}$ (`FP`/`FN`, F1, $\text{MAE}_t$) | AWGN $5\text{ dB}$ (`FP`/`FN`, F1, $\text{MAE}_t$) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`100`** | `0/0`, 100.00%, 1.50 ms | `0/0`, 100.00%, 1.56 ms | `0/0`, 100.00%, 1.53 ms | `0/0`, 100.00%, 1.57 ms | `0/0`, 100.00%, 1.74 ms | `0/0`, 100.00%, 1.96 ms |
| **`101`** | `4/1`, 99.87%, 1.86 ms | `4/1`, 99.87%, 1.84 ms | `4/1`, 99.87%, 1.83 ms | `4/1`, 99.87%, 1.82 ms | `4/1`, 99.87%, 1.92 ms | `4/1`, 99.87%, 2.14 ms |
| **`103`** | `0/1`, 99.98%, 1.92 ms | `0/1`, 99.98%, 1.97 ms | `0/1`, 99.98%, 1.92 ms | `0/1`, 99.98%, 1.96 ms | `0/1`, 99.98%, 2.09 ms | `0/1`, 99.98%, 2.35 ms |
| **`105`** | `33/8`, 99.20%, 1.23 ms | `33/8`, 99.20%, 1.18 ms | `33/7`, 99.22%, 1.65 ms | `36/11`, 99.09%, 1.76 ms | `34/10`, 99.15%, 2.69 ms | `34/12`, 99.11%, **3.56 ms** |
| **`106`** | `0/2`, 99.95%, 2.36 ms | `0/2`, 99.95%, 2.38 ms | `0/2`, 99.95%, 2.80 ms | `0/2`, 99.95%, 2.83 ms | `0/4`, 99.90%, 3.58 ms | `4/6`, 99.75%, **4.27 ms** |

---

## 5. Engineering Answers to Step 16 Diagnostic Questions

1. **Which interference produces the most False Positives (FPs)?**
   - Among *synthetic* perturbations, **severe broadband Gaussian white noise ($5\text{ dB}$ SNR)** and **in-band clinical electrode motion artifacts** (native to Record `105`, producing `33–36 FP`) are the primary drivers of False Positives.
   - Stationary out-of-band interferences ($0.25\text{ Hz}$ respiratory baseline wander and $60\text{ Hz}$ power-line hum) generated **zero additional FPs** on clean and arrhythmic records because their frequencies lie more than two octaves outside the $5\text{–}20\text{ Hz}$ 4th-order zero-phase Butterworth passband (`sosfiltfilt` yields an effective 8th-order $-48\text{ dB/octave}$ stopband attenuation).
   - By contrast, white noise distributes spectral energy uniformly across $[0, 180\text{ Hz}]$; the $\sim 8.3\%$ spectral fraction inside $5\text{–}20\text{ Hz}$ passes unattenuated into the differentiator $d[n]$ and squarer $e[n]$, creating transient energy bursts that cross `THRESHOLD_I1` (e.g., `0 -> 4 FP` on Record `106` at $5\text{ dB}$).

2. **Under which condition do missed beats (False Negatives, FNs) increase?**
   - FNs increase when **high broadband noise ($10\text{ dB}$ and $5\text{ dB}$ SNR) coincides with ectopic or aberrant QRS morphologies** (Record `106`: `2 FN` at baseline $\rightarrow$ `4 FN` at $10\text{ dB}$ $\rightarrow$ `6 FN` at $5\text{ dB}$; Record `105`: `8 FN` $\rightarrow$ `12 FN` at $5\text{ dB}$).
   - **Biophysical mechanism:** Multiform premature ventricular contractions (`V` beats, comprising $25.6\%$ of Record `106`) have wider durations and gentler slopes ($|dV/dt|$) than narrow His–Purkinje-conducted `N` beats. Because the differentiator acts as a high-pass slope amplifier, wide PVCs generate smaller peaks in $e[n] = (dx/dn)^2$ while broadband noise elevates the running noise floor `NPKI`, narrowing the margin between weak ectopic beats and `THRESHOLD_I1`.

3. **How much does Sensitivity change as SNR decreases?**
   - On high-SNR narrow-QRS sinus recordings (`100`, `101`, `103`), Sensitivity remains invariant ($99.95\%\text{–}100.00\%$) down to $5\text{ dB}$ SNR.
   - On records with severe artifacts or multiform ventricular ectopy (`105`, `106`), Sensitivity drops progressively from $99.69\% \rightarrow 99.53\%$ (`105`) and $99.90\% \rightarrow 99.70\%$ (`106`). While aggregate percentages look high, in a 30-minute window those drops represent $50\%\text{–}200\%$ relative increases in missed cardiac beats.

4. **Do detections shift in time ($\text{MAE}_t$ and signed bias)?**
   - **Yes — temporal jitter is a far more sensitive early indicator of signal degradation than beat-count F1-score.**
   - Even on Record `100` where F1 remains $100.00\%$ across all six scenarios, $\text{MAE}_t$ increases monotonically by **$+30.0\%$** (`1.50 ms` $\rightarrow$ `1.96 ms`) and timing error standard deviation increases by **$+56.1\%$** (`1.68 ms` $\rightarrow$ `2.62 ms`).
   - On Record `105` and Record `106`, $\text{MAE}_t$ at $5\text{ dB}$ reaches **`3.56 ms` ($+190\%$)** and **`4.27 ms` ($+81\%$)**, exceeding the single-sample quantization interval $\Delta t = 2.78\text{ ms}$ — a critical degradation if downstream heart-rate variability (HRV) fiducial intervals are computed from detected R-peaks.

5. **Does the $5\text{–}20\text{ Hz}$ filter suppress interference without excessively degrading QRS localization?**
   - Because `scipy.signal.sosfiltfilt` applies forward-backward zero-phase filtering, it eliminates nonlinear phase delay (group delay $= 0$ across all frequencies), maintaining sub-sample localization ($\text{MAE}_t \approx 1.2\text{–}2.4\text{ ms} < \Delta t$) on baseline recordings while completely rejecting $0.25\text{ Hz}$ drift and $60\text{ Hz}$ interference.
   - However, zero-phase `sosfiltfilt` is **non-causal** (requires the full signal buffer or block look-ahead) and cannot run sample-by-sample in real-time embedded patient monitors without introducing causal filter delay.

6. **Does the detector behave identically across distinct recordings?**
   - **No.** Evaluating exclusively on Record `100` (`F1 = 100.0%`, `FP = 0`, `FN = 0`) creates a strong optimistic bias. Expanding to Record `105` (motion artifacts) and Record `106` (multiform PVCs) immediately exposes classical Pan–Tompkins vulnerabilities to **in-band non-stationary electrode motion** ($5\text{–}15\text{ Hz}$ spectral overlap with QRS) and **slope-attenuated ventricular ectopy**.

---

## 6. Known Engineering, Clinical & Metrological Limitations

1. **Not a Medical Device or Diagnostic System:** `ECG-Integrity-Lab` is an offline biomedical engineering verification test bench. It has not undergone clinical validation, risk management (ISO 14971), software lifecycle verification (IEC 62304), or regulatory review.
2. **Non-Diagnostic Filter Bandwidth ($5\text{–}20\text{ Hz}$):** Diagnostic Resting/Ambulatory ECG standards require much wider passbands (typically $0.05\text{–}150\text{ Hz}$ or $0.67\text{–}40\text{ Hz}$ depending on application) to preserve ST-segment elevation/depression, Q-waves, and P/T morphology. The $5\text{–}20\text{ Hz}$ filter intentionally distorts waveform morphology to maximize QRS detection SNR.
3. **Non-Causal Zero-Phase Filtering:** `sosfiltfilt` requires offline batch data; transitioning to a real-time embedded microcontroller requires causal IIR/FIR filtering and explicit group-delay compensation.
4. **Single-Lead Evaluation:** Only Channel 0 (`MLII`) was processed in this baseline study. Multi-lead sensor fusion (combining `MLII` + precordial `V1`/`V5`) is essential when a ventricular ectopic axis is nearly orthogonal to `MLII`.
5. **Synthetic vs. Physiological Noise Models:** Stationary sinusoidal $0.25\text{ Hz}$ drift, stationary $60\text{ Hz}$ hum, and stationary Gaussian white noise do not fully capture non-stationary EMG bursts, abrupt electrode-skin impedance steps (half-cell potential pops), or electrosurgical interference. MIT-BIH Noise Stress Test Database (`nstdb`) electrode-motion (`em`), muscle-artifact (`ma`), and baseline-wander (`bw`) recordings should be incorporated in Phase 2.
6. **Software Phantom vs. Traceable Calibration:** `digital_phantom.py` verifies floating-point software scaling and interval math; it does not verify physical ADC linearity, input impedance, CMRR, or SI-traceable voltage standards.
7. **Dataset Sampling & Statistical Independence:** The 5 pilot records are a purposive engineering stress test, not a random sample. Across the full 48-record MIT-BIH database (47 subjects, since `201` and `202` belong to the same patient), future confidence intervals must use subject-clustered bootstrap resampling.
