# Experimental Protocol & Metrological QA Specification

**Project:** ECG-Integrity-Lab (`v0.1.0`)  
**Subtitle:** Open ECG Signal Conditioning & QRS Detector Verification Framework  
**Document Type:** Pre-Specified Verification & Experimental Protocol  

---

## 1. Research Objective & Technical Question

The primary engineering question investigated in this protocol is:

> **How do controlled electrical interferences and broadband physiological/instrumentation noise degrade the beat-detection accuracy (Sensitivity, PPV, F1-score) and temporal localization precision ($\text{MAE}_t$) of a classical Pan–Tompkins-inspired QRS detection pipeline?**

Rather than claiming clinical diagnostic validity from visual inspection, this protocol establishes a reproducible computer-based test bench combining **biomedical instrumentation principles**, **digital signal conditioning**, **metrological verification**, and **quantitative quality assurance (QA/QC)**.

---

## 2. Dataset & Attribution

- **Source Database:** MIT-BIH Arrhythmia Database (`mitdb`, v1.0.0), hosted on PhysioNet.
- **Database Characteristics:**
  - 48 two-channel ambulatory ECG recordings (~30 minutes each) obtained from **47 subjects** (records `201` and `202` originate from the same male subject; therefore, statistical independence across all 48 records must **not** be assumed).
  - Sampling frequency: $f_s = 360\text{ Hz}$ ($\Delta t = \frac{1}{360}\text{ s} \approx 2.7778\text{ ms}$).
  - Digitization resolution: 11-bit ADC over a $\pm 5\text{ mV}$ ($10\text{ mV}$ peak-to-peak) nominal range (`200` ADC counts/mV, digital baseline `1024` counts in standard records).
  - Expert-reviewed reference annotations (`.atr`): ~110,000 beat annotations.
- **License:** Open Data Commons Attribution License v1.0 (ODC-By v1.0). Full bibliographic attribution to Moody & Mark (2001) and PhysioNet is preserved in [`REFERENCES.md`](REFERENCES.md).

### Staged Record Selection Plan

| Stage | Scope | Selected Records | Engineering Rationale |
| :--- | :--- | :--- | :--- |
| **Level 1 — Initial Verification** | 1 record, first 300 s (5 min) | `100` (Channel 0, `MLII`) | End-to-end verification of WFDB retrieval, header QA, filtering, QRS localization, and controlled noise injection on clean normal sinus rhythm with occasional premature beats. |
| **Level 2 — Multi-Record Pilot** | 5 full records (~1800 s / 30 min each) | `100`, `101`, `103`, `105`, `106` (Channel 0) | Practical engineering stress test across heterogeneous clinical/technical conditions: clean sinus (`100`, `103`), baseline drift/artifact (`101`), severe motion/noise artifacts (`105`), and frequent multiform premature ventricular contractions (`106`). *Note: pragmatic stress-test selection, not a random sample.* |
| **Level 3 — Full Database Study** | All 48 records | `100`–`234` | Population-wide variability analysis accounting for subject clustering (47 subjects) and paced rhythms (`102`, `104`, `107`, `217`). |

---

## 3. Three-Tier Verification & Quality Assurance Architecture

### Tier A — Software & Scale Verification (`digital_phantom.py`)

Before analyzing physiological signals, the mathematical scaling and interval routines are verified against a deterministic digital calibration phantom:

- **Nominal Vector:** Rectangular pulse train of amplitude $A_{\text{ref}} = 1.0\text{ mV}$, pulse width $100\text{ ms}$, rate $60.0\text{ bpm}$ ($T_{RR} = 1.0\text{ s}$), duration $10.0\text{ s}$, sampled at $f_s = 360\text{ Hz}$.
- **Acceptance Criteria:**
  - Detected pulse count $= 10$
  - Amplitude error $E_A = \frac{A_{\text{measured}} - A_{\text{ref}}}{A_{\text{ref}}} \times 100\% = 0.0\%$
  - Rate error $E_R = \frac{\text{BPM}_{\text{measured}} - \text{BPM}_{\text{ref}}}{\text{BPM}_{\text{ref}}} \times 100\% = 0.0\%$
  - Perturbation sensitivity check: an injected $+5\%$ gain factor ($1.05\text{ mV}$) must be reported as $E_A = +5.0\%$.
- **Negative / Fault-Injection Checks (`pytest`):** Explicit rejection (`ValueError`) of invalid sampling frequencies ($f_s \le 0$ or header mismatch), nonexistent channels, `NaN`/`Inf` samples, and uncalibrated raw ADC counts exceeding physiological millivolt limits ($|x| > 20\text{ mV}$).
- **Metrological Disclaimer:** This digital test vector verifies software arithmetic only. It is **not** traceable to the International System of Units (SI), does not replace a physical calibrated ECG simulator, and does not verify hardware analog front-end gain.

### Tier B — Data Integrity & Signal Processing Verification (`dataset.py`, `signal_processing.py`, `detector.py`)

1. **Header & Calibration Audit (Step 9):**
   - Confirm $f_s = 360\text{ Hz}$ and physical unit `mV`.
   - Record header `adc_gain` (counts/mV), digital `baseline` (counts), and `adc_res` (bits). Never process raw ADC integers directly as millivolts.
   - Filter `.atr` annotations strictly against standard PhysioNet QRS beat symbols (`N, L, R, B, A, a, J, S, V, r, F, e, j, n, E, /, f, Q, ?`), explicitly excluding rhythm labels (`+`), signal-quality markers (`~`), isolated artifacts (`|`), and non-conducted P-waves (`x`).
2. **Signal Conditioning Filter (Step 11):**
   - 4th-order Butterworth bandpass ($5\text{–}20\text{ Hz}$) implemented as cascaded second-order sections (`scipy.signal.butter(..., output='sos')`) and applied via zero-phase forward-backward filtering (`scipy.signal.sosfiltfilt`) for numerical stability and zero phase lag.
   - **Bandwidth Scope Warning:** A $5\text{–}20\text{ Hz}$ passband maximizes QRS complex SNR while suppressing respiratory baseline drift and $60\text{ Hz}$ mains pickup, but it **attenuates P/T waves and distorts ST-segment morphology**. Outputs must never be used for diagnostic waveform or ischemia/ST interpretation.
3. **Pan–Tompkins-Inspired QRS Detector Parameters (Step 12):**
   - Numerical derivative: centered gradient $d[n]$.
   - Nonlinear squaring: $e[n] = d[n]^2$.
   - Moving-window integration: $M[n] = \frac{1}{N}\sum_{k=0}^{N-1} e[n-k]$ over a $120\text{ ms}$ window ($N = 43$ samples at $360\text{ Hz}$).
   - Adaptive dual-threshold peak detection (`SPKI`, `NPKI`, `THRESHOLD_I1`, searchback `THRESHOLD_I2` at $1.66 \times \overline{RR}$) with $250\text{ ms}$ physiological refractory period ($90$ samples).
   - QRS peak relocation: local extremum search within $\pm 90\text{ ms}$ ($\pm 32$ samples) enforcing $250\text{ ms}$ minimum inter-beat separation.

### Tier C — Controlled Perturbations & Quantitative Evaluation (`evaluation.py`)

Each ECG segment is evaluated across six controlled experimental conditions:

| # | Scenario Key | Interference Type | Experimental Specification |
| :-: | :--- | :--- | :--- |
| 1 | `baseline` | Unperturbed reference | Calibrated physical ECG signal $x[n]$ |
| 2 | `baseline_wander` | Respiratory / motion drift | Sinusoidal $0.20\text{ mV}$ peak at $0.25\text{ Hz}$ |
| 3 | `mains_interference` | Power-line hum | Sinusoidal $0.05\text{ mV}$ peak at $60.0\text{ Hz}$ |
| 4 | `white_noise_20dB` | Broadband Gaussian noise | Zero-mean AWGN at $\text{SNR} = 20\text{ dB}$ (`seed=42`) |
| 5 | `white_noise_10dB` | Broadband Gaussian noise | Zero-mean AWGN at $\text{SNR} = 10\text{ dB}$ (`seed=42`) |
| 6 | `white_noise_5dB` | Broadband Gaussian noise | Zero-mean AWGN at $\text{SNR} = 5\text{ dB}$ (`seed=42`) |

Signal-to-noise ratio is referenced to AC signal variance $P_s = \text{Var}(x)$ so static DC offsets cannot bias injected noise amplitude:
$$P_n = \frac{\text{Var}(x)}{10^{\text{SNR}_{\text{dB}} / 10}}, \qquad \sigma_n = \sqrt{P_n}$$

#### One-to-One Beat Matching & Performance Metrics

- **Matching Window:** $\pm 150\text{ ms}$ ($\pm 54$ samples at $360\text{ Hz}$), using greedy minimum-distance **strict one-to-one assignment** so multiple detections near a single beat cannot inflate True Positives.
- **Classification & Timing Metrics:**
  $$\text{Sensitivity (Se)} = \frac{TP}{TP + FN}, \qquad \text{PPV} = \frac{TP}{TP + FP}, \qquad F_1 = \frac{2 \cdot \text{Se} \cdot \text{PPV}}{\text{Se} + \text{PPV}}$$
  $$\text{MAE}_t = \frac{1}{N_{TP}} \sum_{i=1}^{N_{TP}} |t_{\text{det}, i} - t_{\text{ref}, i}| \quad [\text{ms}]$$

---

## 4. Regulatory & Standards Context

- **ANSI/AAMI EC57:2012 (R2020)** (*Testing and reporting performance results of cardiac rhythm and ST segment measurement algorithms*) provides industry context for beat-by-beat evaluation windows and performance reporting. This open framework implements an educational/research beat-matching evaluator inspired by those principles; it is **not** a complete implementation of ANSI/AAMI EC57 (e.g., it omits the standard's 5-minute learning-period exclusion rules, ventricular flutter episode handling, and multi-database mandatory test sets) and **claims no regulatory conformity**.
- **IEC 60601-2-47:2012** addresses basic safety and essential performance of ambulatory electrocardiographic systems (hardware + embedded firmware). Because `ECG-Integrity-Lab` is a purely offline software pipeline without patient-connected hardware, electrical safety and medical-device essential performance conformity are outside its scope.
