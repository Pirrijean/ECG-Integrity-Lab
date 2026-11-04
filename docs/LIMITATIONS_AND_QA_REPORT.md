# Quantitative QA Report, Phase-2 Error Analysis & Methodological Limitations (`v0.2.0`)

**Project:** ECG-Integrity-Lab (`v0.2.0`)  

---

## 1. Pooled Baseline Performance (Micro-Aggregated Counts vs. Per-Record Breakdown)

Across the five 30-minute MIT-BIH pilot recordings (`100`, `101`, `103`, `105`, `106`), pooling counts directly (micro-aggregation, avoiding unweighted averaging of percentages) yields:

| Metric | Pooled Baseline Count / Value | Methodological Note |
| :--- | :---: | :--- |
| **Reference Beats ($N_{\text{ref}}$)** | `10,789` | True QRS beat symbols only (`190+` non-beat `.atr` markers excluded) |
| **True Positives (TP)** | `10,777` | Strict 1-to-1 matching ($\pm 150\text{ ms}$ primary criterion) |
| **False Positives (FP)** | `37` | **89.2% (`33/37`) concentrated in Record `105`** |
| **False Negatives (FN)** | `12` | **`8/12` in Record `105`; `2/12` in Record `106`** |
| **Pooled Sensitivity** | **99.89%** | $10,777 / (10,777 + 12)$ |
| **Pooled PPV** | **99.66%** | $10,777 / (10,777 + 37)$ |
| **Pooled F1-Score** | **99.77%** | Harmonic mean of pooled Se and PPV |

> [!NOTE]
> Evaluating these `10,789` beats across `6` perturbation conditions produces `64,734` *repeated-measures evaluation opportunities*, **not** `64,734` statistically independent observations.

---

## 2. Core Finding: F1-Score Stability vs. Temporal Jitter & HRV (`RMSSD`) Bias

On Record `100` (30 min, $2,265$ reference beats; discrete sampling interval $\Delta t = 2.78\text{ ms}$), `TP = 2,265`, `FP = 0`, `FN = 0` (`F1 = 100.00%`) across all synthetic SNR levels, yet temporal precision and downstream normal-to-normal (`N`–`N`, $n = 2,196$ clean sinus pairs) `RMSSD` degrade monotonically (`results/pilot/record_100_ch0/hrv_jitter_metrics.csv`):

| Condition | F1-Score (%) | Ensemble $\text{MAE}_t$ (ms) | Error SD $\sigma_t$ (ms) | RR MAE (ms) | Reference `RMSSD` (ms) | Detected `RMSSD` (ms) | `RMSSD` Inflation (%) |
| :--- | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
| **Baseline** | 100.00 | 1.50 | 1.68 | 1.52 | 27.80 | 28.10 | +1.06% |
| **AWGN 20 dB** | 100.00 | 1.57 | 1.73 | 1.57 | 27.80 | 28.06 | +0.94% |
| **AWGN 10 dB** | 100.00 | 1.74 | 2.31 | 1.93 | 27.80 | 28.27 | +1.68% |
| **AWGN 5 dB** | 100.00 | **1.96** (+30.7%) | **2.62** (+56.0%) | **2.41** | 27.80 | **28.51** | **+2.54%** |

*(Precision note: $\text{MAE}_t < 2.78\text{ ms}$ reflects the arithmetic average of discrete sample offsets across $2,265$ beats; individual peak indices remain quantized to multiples of $\Delta t = 2.78\text{ ms}$. Signed errors $e_i = t_{\text{det},i} - t_{\text{ref},i} > 0$ indicate detections slightly lag the expert `.atr` mark).*

---

## 3. Experiment 1 — Sensitivity to Temporal Matching Window (`tolerance_sweep.csv`)

Evaluating Record `106` ($2,021$ beats, $25.6\%$ multiform PVCs) across six matching tolerances demonstrates why a single $\pm 150\text{ ms}$ window overstates fiducial accuracy:

| Matching Window | Samples ($360\text{ Hz}$) | Baseline F1 (%) (`FP`/`FN`) | AWGN $10\text{ dB}$ F1 (%) (`FP`/`FN`) | AWGN $5\text{ dB}$ F1 (%) (`FP`/`FN`) |
| :---: | :---: | :---: | :---: | :---: |
| **`±8.33 ms`** | `3` | **95.79%** (`84/86`) | **95.74%** (`84/88`) | **93.96%** (`121/123`) |
| **`±16.67 ms`** | `6` | 97.13% (`57/59`) | 97.33% (`52/56`) | 97.03% (`59/61`) |
| **`±25.0 ms`** | `9` | 97.38% (`52/54`) | 97.37% (`51/55`) | 97.23% (`55/57`) |
| **`±50.0 ms`** | `18` | 99.80% (`3/5`) | 99.80% (`2/6`) | 99.60% (`7/9`) |
| **`±100.0 ms`** | `36` | 99.95% (`0/2`) | 99.90% (`0/4`) | 99.70% (`5/7`) |
| **`±150.0 ms`** | `54` | 99.95% (`0/2`) | 99.90% (`0/4`) | 99.75% (`4/6`) |

Approximately `50–85` multiform PVCs are detected `25–50 ms` away from the expert annotation (where the local filtered extremum lands on a broad S-wave or secondary R'-lobe rather than the clinician's fiducial mark), causing a sharp drop in F1 below $\pm 50\text{ ms}$.

---

## 4. Experiment 2 — Real Physiological Acquisition Noise (`nstdb`: `bw`, `ma`, `em`)

Testing Record `100` (30 min, $2,265$ beats) against physical noise recordings from the **MIT-BIH Noise Stress Test Database** (`results/nstdb/record_100_ch0/nstdb_metrics.csv`) highlights the gap between synthetic sinusoids/AWGN and ambulatory artifacts:

| Physical Noise Source (`nstdb`) | $24\text{ dB}$ F1 | $12\text{ dB}$ F1 | $6\text{ dB}$ F1 | $0\text{ dB}$ F1 (`FP`/`FN`) | $-6\text{ dB}$ F1 (`FP`/`FN`, PPV, $\sigma_t$) |
| :--- | :-: | :-: | :-: | :-: | :-: |
| **`bw` (Baseline Wander)** | 100.00% | 100.00% | 100.00% | 100.00% (`0/0`) | **100.00%** (`0/0`, PPV 100.0%, 1.69 ms) |
| **`ma` (Muscle / EMG Artifact)** | 100.00% | 100.00% | 100.00% | 99.98% (`1/0`) | **97.18%** (`127/4`, PPV **94.68%**, 5.25 ms) |
| **`em` (Electrode Motion)** | 100.00% | 100.00% | 100.00% | 100.00% (`0/0`) | **96.27%** (`171/4`, PPV **92.97%**, **6.38 ms**) |

Furthermore, adding `nstdb/em` at $+6\text{ dB}$ onto Record `105` increased False Positives from `33` to `43` (`PPV = 98.34%`). Physical electrode motion (`em`) is substantially more destructive to PPV than stationary Gaussian noise because mechanical skin–electrode transients directly mimic QRS complexes inside $5\text{–}20\text{ Hz}$.

---

## 5. Experiment 3 — Empirical Verification of Record `105` FPs & Record `106` PVC FNs

### 5.1 Record `105` False Positive Inspection (`results/diagnostics/record_105_baseline/fp_inspection.csv`)
Inspecting $\pm 300\text{ ms}$ windows around all `33` Baseline False Positives confirmed:
- **Direct Artifact Co-Localization:** FPs `#1`, `#2`, and `#3` (`t = 631.45 s, 1029.94 s, 1030.79 s`) occur within **`3 ms` to `36 ms`** of expert PhysioNet `|` (*isolated QRS-like artifact*) annotations.
- **Severe Signal-Quality Bursts:** `19/33` FPs occur within $5\text{ s}$ of PhysioNet `~` (*signal-quality/noise transition*) markers, and the remaining `11` lie inside the prolonged `1213–1245 s` high-amplitude motion episode (`raw_ptp_mv` between `2.35` and `5.77 mV`), where $5\text{–}20\text{ Hz}$ energy surges to **`1.57×`–`27.74×` `THRESHOLD_I1`**.

### 5.2 Record `106` Biophysical PVC Comparison (`results/diagnostics/record_106_awgn_5dB/`)
Stratifying Record `106` at $5\text{ dB}$ SNR by beat symbol showed **`100.00%` Sensitivity on normal `N` beats (`1503/1503`, `0 FN`)** vs. **`98.84%` Sensitivity on `V` beats (`512 TP`, `6 FN`)**, with PVC $\text{MAE}_t = 8.28\text{ ms}$ vs. $2.91\text{ ms}$ on `N` beats. Comparing the four biophysical parameters across Detected (`512 TP`) vs. Missed (`6 FN`) PVCs empirically confirmed the energy-attenuation mechanism:
- **Peak-to-peak amplitude $A_{\text{QRS}}$:** `TP` mean $\approx 2.57\text{ mV}$ vs. `FN` mean $\approx 1.50\text{ mV}$.
- **Filtered slope $\max |dV/dt|$:** `TP` mean $\approx 58.8\text{ mV/s}$ vs. **`FN` mean $\approx 18.6\text{ mV/s}$ ($3.2\times$ lower slope)**.
- **Normalized integrated energy $E_{\text{QRS}} / I_1$:** `TP` mean $= 2.89\times I_1$, whereas **all `6/6` missed PVCs had $E_{\text{QRS}} / I_1 < 0.60$** (strictly below the adaptive detection threshold).

---

## 6. Experiment 4 — Independent Benchmark vs. WFDB `XQRS` (`results/benchmark/`)

Comparing `Pan–Tompkins (Custom)` against `wfdb.processing.xqrs_detect` under identical inputs ($\pm 150\text{ ms}$ primary window) revealed a clear engineering trade-off rather than uniform superiority of either method:

| Record (30 min) | Condition | Custom Pan–Tompkins (`FP`/`FN`, F1, $\text{MAE}_t$) | WFDB `XQRS` (`FP`/`FN`, F1, $\text{MAE}_t$) | Key Engineering Insight |
| :---: | :--- | :---: | :---: | :--- |
| **`100`** | Baseline | `0/0`, 100.00%, 1.50 ms | `0/0`, 100.00%, **0.50 ms** | XQRS achieves tighter R-peak centering on narrow sinus beats |
| **`105`** | Baseline | `33/8`, 99.20%, 1.23 ms | `33/4`, **99.28%**, **1.16 ms** | Both suffer `33 FP` on in-band motion; XQRS halves `FN` (`8 -> 4`) |
| **`105`** | Real EM $6\text{ dB}$ | `43/10`, 98.97%, 1.45 ms | `34/4`, **99.26%**, **1.27 ms** | XQRS Ricker wavelet rejects `+9` electrode-motion FPs |
| **`106`** | Baseline | **`0/2`, 99.95%**, 2.36 ms | `0/88`, 97.77%, 1.25 ms | **XQRS misses `88` multiform PVCs (`Se = 95.65%`) vs. `2` for Custom PT** |
| **`106`** | AWGN $5\text{ dB}$ | **`4/6`, 99.75%**, 4.27 ms | `22/93`, 97.10%, 2.13 ms | Custom PT retains higher PVC sensitivity; XQRS keeps lower jitter on TPs |

---

## 7. Summary of Limitations

1. **Exploratory Characterization, Not Clinical Validation:** Results come from 5 purposively chosen MIT-BIH records; full 48-record evaluation with subject-clustered bootstrap CIs (47 subjects) is required before generalizing across ambulatory populations.
2. **Non-Diagnostic $5\text{–}20\text{ Hz}$ Bandwidth & Non-Causal `sosfiltfilt`:** Unsuitable for ST-segment or P/T morphological diagnosis and requires causal redesign with group-delay compensation for real-time embedded hardware.
3. **Functional Digital Phantom Only:** Verifies software scale/rate math (`0.0%` nominal error, `+5.0%` at $1.05\times$), not hardware analog front-ends or SI traceability.
