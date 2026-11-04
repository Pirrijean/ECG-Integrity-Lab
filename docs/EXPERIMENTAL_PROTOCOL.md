# Experimental Protocol & Verification Specification (`v0.2.0`)

**Project:** ECG-Integrity-Lab (`v0.2.0`)  
**Subtitle:** Open ECG Signal Conditioning & QRS Detector Verification Framework  
**Document Type:** Pre-Specified Experimental & Methodological Protocol  

---

## 1. Research Objective & Central Thesis

The central engineering thesis investigated by this framework is:

> **Beyond F1-Score: Beat-detection accuracy (Sensitivity, PPV, F1) can remain near-perfect while R-peak temporal localization precision ($\text{MAE}_t$, $\sigma_t$) and downstream time-domain HRV metrics (`RMSSD`) already undergo measurable degradation—and detection failure modes depend strongly on whether disturbances are synthetic out-of-band sinusoids, broadband noise, physical electrode motion, or multiform ventricular ectopy.**

Rather than claiming that a classical detector is universally "robust", `ECG-Integrity-Lab` provides an open, reproducible computer-based bench to **characterize where, how, and why** an ECG signal-processing pipeline degrades.

---

## 2. Datasets & Attribution

1. **MIT-BIH Arrhythmia Database (`mitdb`, v1.0.0 — PhysioNet):**
   - 48 two-channel ambulatory recordings (~30 min each) from **47 subjects** (`201` and `202` belong to the same patient; records are therefore **not** statistically independent across all 48 files).
   - $f_s = 360\text{ Hz}$ ($\Delta t = \frac{1}{360}\text{ s} \approx 2.7778\text{ ms}$), 11-bit ADC over $\pm 5\text{ mV}$ (`200 counts/mV`, digital baseline `1024`).
   - Five-record purposive stress-test pilot (`100`, `101`, `103`, `105`, `106`; `10,789` true reference beats): selected deliberately to contrast clean sinus rhythm (`100`, `103`), mild drift/artifact (`101`), severe motion/artifact (`105`), and frequent multiform PVCs (`106`). *Not a random sample of MIT-BIH.*
2. **MIT-BIH Noise Stress Test Database (`nstdb`, v1.0.0 — PhysioNet):**
   - Physical ambulatory noise recordings ($f_s = 360\text{ Hz}$, `mV`) acquired from active volunteers with electrodes placed on the limbs/thorax where ECG is absent:
     - `bw`: physical baseline wander
     - `ma`: physical muscle (EMG) artifact
     - `em`: physical electrode motion artifact

---

## 3. Verification Architecture & Methodological Conventions

### 3.1 Functional and Mathematical Verification via Digital Reference Signals (`digital_phantom.py`)
- Uses a deterministic synthetic rectangular pulse train ($1.0\text{ mV}$, $60.0\text{ bpm}$, $100\text{ ms}$ width, $10\text{ s}$ at $360\text{ Hz}$) plus a $+5\%$ gain perturbation ($1.05\text{ mV} \rightarrow E_A = +5.0\%$).
- **Terminology Precision:** This is a **functional and mathematical software verification test**, *not* a metrological hardware calibration. It does not verify an electrocardiograph's analog front-end, ADC linearity, or traceability to SI units.

### 3.2 Temporal Resolution vs. Ensemble $\text{MAE}_t$ & Signed Error Convention (`evaluation.py`)
- **Discrete Sampling Resolution:** Because `detect_qrs_pan_tompkins` outputs integer sample indices at $f_s = 360\text{ Hz}$, the discrete resolution of any single beat detection is $\Delta t = 2.7778\text{ ms}$. Ensemble averages such as $\text{MAE}_t = 1.17\text{ ms}$ or $1.50\text{ ms}$ arise mathematically from averaging discrete sample offsets ($0, \pm 1, \pm 2, \dots$ samples) across hundreds or thousands of beats; they must **never** be interpreted as sub-sample single-beat resolution unless explicit sub-sample interpolation (`refine_peaks_parabolic`) is enabled and independently validated.
- **Signed Timing Error Convention:**
  $$e_i = t_{\text{det}, i} - t_{\text{ref}, i} = \frac{n_{\text{det}, i} - n_{\text{ref}, i}}{f_s}$$
  A **positive** `mean_signed_error_ms` indicates that detections occur systematically *after* the expert `.atr` fiducial mark; a **negative** value indicates detections *precede* the reference.

### 3.3 Sensitivity to Temporal Matching Window (Experiment 1)
Reporting only a permissive $\pm 150\text{ ms}$ ($\pm 54$ samples) window can mask large fiducial shifts. Every evaluation therefore computes strict **one-to-one beat matching** across six tolerances:

| Tolerance ($\pm\text{ms}$) | Samples at $360\text{ Hz}$ | Interpretation |
| :---: | :---: | :--- |
| `±8.33 ms` | `3` | Strict R-peak fiducial precision ($\le 3$ samples) |
| `±16.67 ms` | `6` | Tight HRV-grade R-peak alignment |
| `±25.0 ms` | `9` | Narrow QRS half-width window |
| `±50.0 ms` | `18` | Moderate QRS complex width window |
| `±100.0 ms` | `36` | Wide ectopic QRS tolerance |
| `±150.0 ms` | `54` | Primary beat-presence tolerance (reference benchmark) |

### 3.4 Real Physiological Noise Mixing Protocol (Experiment 2 — `nstdb.py`)
Clean ECG $x[n]$ and zero-mean physical noise $n_0[n] = n[n] - \bar{n}$ are combined at target $\text{SNR} \in \{24, 18, 12, 6, 0, -6\}\text{ dB}$ while preserving original `.atr` beat timestamps:
$$x_{\text{noisy}}[n] = x[n] + \alpha \, n_0[n], \qquad \alpha = \sqrt{\frac{P_x}{P_n \cdot 10^{\text{SNR}_{\text{dB}} / 10}}}$$
where $P_x = \text{Var}(x)$ and $P_n = \text{Var}(n_0)$ are computed over the exact aligned recording segment.

### 3.5 Morphology-Stratified & Biophysical Error Diagnostics (Experiment 3 — `error_analysis.py`)
- Stratifies Sensitivity and $\text{MAE}_t$ by PhysioNet beat symbol (`N`, `A`, `V`, `Q`).
- Extracts four biophysical/signal parameters around every reference beat ($\pm 80\text{ ms}$): peak-to-peak amplitude $A_{\text{QRS}}$, filtered half-maximum width $W_{\text{QRS}}$, peak slope $\max |dV/dt|$, and moving-window integrated energy $E_{\text{QRS}}$ relative to the local adaptive threshold `THRESHOLD_I1`.
- Extracts $\pm 300\text{ ms}$ multi-stage windows (raw, $5\text{–}20\text{ Hz}$ filtered, derivative $d[n]$, integrated $M[n]$ + `THRESHOLD_I1`) around False Positives and cross-references auxiliary `.atr` markers (`~`, `|`, `+`).

### 3.6 Independent Detector Benchmark (Experiment 4 — `benchmark.py`)
Compares `detect_qrs_pan_tompkins` against `wfdb.processing.xqrs_detect` (`learn=True`) on identical signals, reference annotations, and tolerance sweeps. *(Note: both operate on the $5\text{–}20\text{ Hz}$ QRS band, isolating differences between derivative-squared moving average + dual Pan–Tompkins thresholds vs. XQRS Ricker-wavelet envelope integration and learning thresholds).*

### 3.7 Time-Domain HRV Localization-Jitter Propagation (`hrv.py`)
Isolates verified consecutive normal sinus (`N`–`N`) TP beat pairs ($300\text{ ms} \le RR_i \le 2000\text{ ms}$), excluding ectopic (`A`, `V`, `Q`) and missed/spurious transitions, to quantify how timing errors $RR_i' = RR_i + e_{i+1} - e_i$ bias `RMSSD` and `SDNN`.

### 3.8 Statistical Structure & Aggregation Rules
- **Micro-Aggregation vs. Macro-Averaging:** Pooled baseline metrics sum `TP`, `FP`, and `FN` across records before computing Sensitivity, PPV, and F1 (`compute_micro_aggregated_metrics`), rather than averaging per-record percentages of unequal beat counts.
- **Repeated Measures Warning:** Evaluating $10,789$ beats under $6$ perturbation conditions yields $64,734$ *evaluation opportunities*, **not** $64,734$ independent observations. Future confidence intervals must use paired differences clustered by record/subject (47 unique subjects across the 48 MIT-BIH records).
