# Experimental Protocol & Verification Specification (`v0.3.0`)

**Project:** ECG-Integrity-Lab (`v0.3.0`)  
**Subtitle:** Open ECG Signal Conditioning, Temporal Metrology & Evidence-Driven QRS Verification Framework  
**Document Type:** Pre-Specified Three-Stage Experimental & Methodological Protocol  

---

## 1. Research Objective & Three-Stage Scientific Progression

The central engineering thesis investigated across the three versions of this framework is:

> **Beyond F1-Score: Beat-detection accuracy (Sensitivity, PPV, F1) under a permissive `±150 ms` window can remain near `100%` while R-peak temporal localization precision (`MAE_t`, `σ_t`) and downstream time-domain HRV metrics (`RMSSD`) already degrade—and once specific failure mechanisms (raw-signal snap jitter, in-band motion bursts, and compensatory-pause PVC masking) are isolated empirically (`v0.2.0`), targeted signal-processing enhancements (`Pan–Tompkins v2` in `v0.3.0`) can be verified reproducibly against both the frozen baseline (`v1`) and an external reference (`WFDB XQRS`).**

### Evolutionary Protocol Stages
- **Stage 1 (`v0.1.0` — Baseline Conditioning & Perturbation Bench):** Verify ADC scaling/header metadata, validate mathematical software routines via a digital reference pulse train, and evaluate classical zero-phase `5–20 Hz` Butterworth + `Pan–Tompkins v1` under controlled synthetic perturbations (`0.25 Hz` wander, `60 Hz` mains, `20/10/5 dB` AWGN).
- **Stage 2 (`v0.2.0` — Critical Failure-Mode & Metrology Characterization):** Freeze `Pan–Tompkins v1` and introduce 6-window tolerance sweeps (`±8.33` to `±150 ms`), physical `NSTDB` noise stress tests (`bw`, `ma`, `em` from `24` to `-6 dB`), `±300 ms` False Positive waveform galleries, biophysical `TP` vs. `FN` feature comparison, `N–N` `RMSSD` jitter propagation, and 2-way benchmarking against `WFDB XQRS`.
- **Stage 3 (`v0.3.0` — Evidence-Driven Algorithmic Upgrade & Three-Way Verification):** Implement `detect_qrs_pan_tompkins_v2` (`v2`) targeting the three exact failure modes diagnosed in Stage 2 while preserving `detect_qrs_pan_tompkins` (`v1`) unchanged, and run a paired 3-way comparison (`PT v1` vs. `PT v2` vs. `WFDB XQRS`).

---

## 2. Datasets & Attribution

1. **MIT-BIH Arrhythmia Database (`mitdb`, v1.0.0 — PhysioNet):**
   - 48 two-channel ambulatory recordings (~30 min each) from **47 unique subjects** (`201` and `202` belong to the same patient; records are therefore **not** statistically independent across all 48 files).
   - Sampling frequency `fs = 360 Hz` (`Δt = 1/360 s ≈ 2.7778 ms`), 11-bit ADC over `±5 mV` (`200 counts/mV`, digital baseline `1024`).
   - Five-record purposive stress-test pilot (`100`, `101`, `103`, `105`, `106`; `10,789` true reference beats): chosen deliberately to contrast clean sinus rhythm (`100`, `103`), mild baseline drift (`101`), severe motion/tape artifact (`105`), and frequent multiform PVCs (`106`, `25.6% V`). *Not a random sample of MIT-BIH.*
2. **MIT-BIH Noise Stress Test Database (`nstdb`, v1.0.0 — PhysioNet):**
   - Physical ambulatory noise recordings (`fs = 360 Hz`, `mV`) acquired from volunteers with electrodes placed on limbs/thorax regions devoid of ECG signal:
     - `bw`: physical baseline wander
     - `ma`: physical muscle (EMG) artifact
     - `em`: physical electrode motion artifact

---

## 3. Verification Architecture & Methodological Conventions

### 3.1 Functional and Mathematical Verification via Digital Reference Signals (`digital_phantom.py`)
- Uses a deterministic synthetic rectangular pulse train (`1.0 mV`, `60.0 bpm`, `100 ms` pulse width, `10 s` at `360 Hz`) plus a `+5%` gain perturbation (`1.05 mV -> E_A = +5.0%`).
- **Terminology Precision:** This is a **functional and mathematical software verification test**, *not* a metrological hardware calibration. It does not calibrate physical ECG front-end amplifiers, ADC linearity, or SI traceability.

### 3.2 Temporal Resolution, Sub-Sample Interpolation & Signed Error Convention (`evaluation.py`)
- **Discrete Sampling Resolution:** Integer sample detections at `fs = 360 Hz` have a single-beat discrete step `Δt = 2.7778 ms`. Ensemble averages below `2.78 ms` (such as `0.27 ms` or `1.50 ms`) arise mathematically from averaging discrete sample offsets (`0, ±1, ±2, ...` samples) across hundreds or thousands of beats.
- **Parabolic Sub-Sample Vertex Interpolation (`v0.3.0`):** `detect_qrs_pan_tompkins_v2` fits a 3-point parabola around `(n0 - 1, n0, n0 + 1)` on `|y_5_20[n]|`:

$$
\delta = \frac{y[n_0 - 1] - y[n_0 + 1]}{2 \left( y[n_0 - 1] - 2 y[n_0] + y[n_0 + 1] \right)}, \qquad n_{\mathrm{sub}} = n_0 + \delta
$$

- **Signed Timing Error Convention:**

$$
e_i = t_{\mathrm{det},i} - t_{\mathrm{ref},i} = \frac{n_{\mathrm{det},i} - n_{\mathrm{ref},i}}{f_s}
$$

A **positive** `mean_signed_error_ms` indicates detections occur systematically *after* the expert `.atr` fiducial mark; a **negative** value indicates detections *precede* the `.atr` mark.

### 3.3 Sensitivity to Temporal Matching Window (Experiment 1)
Every evaluation computes strict **one-to-one beat matching** across six tolerances:

| Tolerance (`±ms`) | Samples at `360 Hz` | Metrological Interpretation |
| :---: | :---: | :--- |
| `±8.33 ms` | `3` | Strict R-peak fiducial precision (`<= 3` samples) |
| `±16.67 ms` | `6` | Tight HRV-grade R-peak alignment |
| `±25.0 ms` | `9` | Narrow QRS half-width window |
| `±50.0 ms` | `18` | Moderate QRS complex width window |
| `±100.0 ms` | `36` | Wide ectopic QRS tolerance |
| `±150.0 ms` | `54` | Primary beat-presence tolerance (reference benchmark) |

### 3.4 Real Physiological Noise Mixing Protocol (Experiment 2 — `nstdb.py`)
Clean ECG `x[n]` and zero-mean physical noise `n0[n] = n[n] - mean(n)` are combined at target `SNR` in `{24, 18, 12, 6, 0, -6} dB` while preserving original `.atr` beat timestamps:

$$
x_{\mathrm{noisy}}[n] = x[n] + \alpha \, n_0[n], \qquad \alpha = \sqrt{\frac{P_x}{P_n \cdot 10^{\mathrm{SNR}_{\mathrm{dB}} / 10}}}
$$

where `P_x = Var(x)` and `P_n = Var(n0)` over the exact aligned segment.

### 3.5 Morphology-Stratified & Biophysical Error Diagnostics (Experiment 3 — `error_analysis.py`)
- Stratifies Sensitivity and `MAE_t` by PhysioNet beat symbol (`N`, `A`, `V`, `Q`).
- Extracts four biophysical parameters around every reference beat (`±80 ms`): peak-to-peak amplitude `A_QRS`, filtered half-maximum width `W_QRS`, peak slope `max|dV/dt|`, and integrated energy ratio `E_QRS / THRESHOLD_I1`.
- Extracts `±300 ms` multi-stage waveforms (raw, `5–20 Hz` filtered, derivative `d[n]`, integrated `M[n]` + `THRESHOLD_I1`) around False Positives and cross-references auxiliary `.atr` markers (`~`, `&#124;`, `+`).

### 3.6 `Pan–Tompkins v2` Design Specification (`v0.3.0` — `detector.py`)
1. **Target 1 (Conditioned Apex Relocation + Parabolic Sub-Sample Vertex):** Eliminates `v1`'s `±15 ms` raw-signal snap so high-frequency noise (`>20 Hz`) cannot perturb `argmax`, relocating directly on the zero-phase `5–20 Hz` extremum and computing `r_peaks_subsample`.
2. **Target 2 (Same-Band `QFSNR` Motion Artifact Gate):** Computes peak amplitude in `±60 ms` divided by flank RMS in `[-240..-85] U [+85..+240] ms` strictly within `5–20 Hz` (`qf5`). Excluding the `±85 ms` QRS core avoids penalizing wide `130 ms` PVCs while rejecting abrupt short-RR (`<0.52 * median_RR`) clutter bursts (`qf5 < 2.8`) and `200–360 ms` post-QRS motion doublets (`qf5 < 3.4` with `<35%` integrated energy of the preceding QRS).
3. **Target 3 (Compensatory-Pause-Aware `2.5–15 Hz` Ectopic Rescue):** Estimates clean `sinus_RR` (`0.75–1.35 * median_RR`), distinguishes true post-PVC compensatory pauses (`prev_RR < 0.72 * sinus_RR`) from genuine `2-RR` missed-beat gaps (`> 1.56 * sinus_RR`), and rescues low-slope PVCs inside genuine `2-RR` gaps using a `2.5–15 Hz` zero-phase envelope gated by same-band `QFSNR >= 3.2`. Optional `secondary_signal` (`Channel 1`) fusion provides backup rescue when dual leads are supplied.

### 3.7 Independent Detector Benchmark & HRV Jitter Propagation (`benchmark.py`, `hrv.py`)
- Compares `Pan–Tompkins v1`, `Pan–Tompkins v2`, and `wfdb.processing.xqrs_detect` (`learn=True`) on identical signals, reference annotations, and tolerance sweeps.
- Isolates verified consecutive normal sinus (`N–N`) TP beat pairs (`300 ms <= RR_i <= 2000 ms`) to quantify how discrete and parabolic sub-sample timing errors `RR'_i = RR_i + e_{i+1} - e_i` bias `RMSSD` and `SDNN`.

### 3.8 Statistical Structure & Aggregation Rules
- **Micro-Aggregation vs. Macro-Averaging:** Pooled summaries sum `TP`, `FP`, and `FN` across records before computing Sensitivity, PPV, and F1 (`compute_micro_aggregated_metrics`).
- **Repeated Measures Warning:** Evaluating `10,789` beats across `6–7` perturbation conditions yields paired repeated-measures evaluation opportunities, **not** independent observations. Inferential extensions across all 48 MIT-BIH records must cluster bootstrap confidence intervals by unique subject (47 subjects).
