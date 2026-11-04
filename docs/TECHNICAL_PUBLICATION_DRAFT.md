# Beyond F1-Score: Evaluating and Improving Temporal Accuracy, Artifact Rejection, and Ectopic Sensitivity in QRS Detection

*A three-stage reproducible investigation using the MIT-BIH Arrhythmia and Noise Stress Test Databases*

**Author:** Pirrijean  
**Repository:** [https://github.com/Pirrijean/ECG-Integrity-Lab](https://github.com/Pirrijean/ECG-Integrity-Lab) (`v0.3.0`)  

---

## Abstract

Evaluations of electrocardiogram (ECG) QRS detectors commonly rely on aggregate beat-classification metrics (Sensitivity, Positive Predictive Value, and F1-score) under a single `±150 ms` matching tolerance and synthetic sinusoidal/Gaussian noise. Using **ECG-Integrity-Lab**, an open Python verification framework, we conduct a three-stage investigation across five 30-minute recordings from the MIT-BIH Arrhythmia Database (`10,789` expert-annotated beats), physical ambulatory noise from the MIT-BIH Noise Stress Test Database (`nstdb`: `bw`, `ma`, `em` from `24` to `-6 dB`), six matching tolerances (`±8.33` to `±150 ms`), downstream normal-to-normal (`N–N`) Heart Rate Variability (`RMSSD`), and `WFDB XQRS`. First (`v0.1.0`–`v0.2.0`), we show that **beat-classification F1 can remain at `100.00%` while temporal R-peak localization degrades substantially**: on `Record 100`, `5 dB` white noise caused zero false positives or missed beats in a classical Pan–Tompkins (`v1`) detector yet increased ensemble timing error (`MAE_t`) by `+30.7%` (`1.50 → 1.96 ms`) and inflated sinus `RMSSD` by `+2.54%`. Biophysical diagnostics revealed three root causes of algorithmic error: (1) unfiltered raw-signal R-peak snapping (`1.50 ms` baseline `MAE_t`), (2) in-band `5–20 Hz` motion bursts (`33` FPs on `Record 105`; `127` `ma` and `171` `em` FPs at `-6 dB`), and (3) `3.2x` lower `dV/dt` slope plus compensatory-pause masking on multiform PVCs (`Record 106`, where `WFDB XQRS` missed `88` beats). Second (`v0.3.0`), we introduce **`Pan–Tompkins v2`**, combining zero-phase `5–20 Hz` apex relocation with parabolic sub-sample interpolation, a core-excluded QRS-to-Flank SNR (`QFSNR`) artifact gate, and compensatory-pause-aware `2.5–15 Hz` ectopic rescue. Across all `10,789` pooled beats, `Pan–Tompkins v2` reduced total classification errors (`FP + FN`) by **`-26.5%` at baseline (`36` vs. `49` for `v1` and `132` for `XQRS`)** and **`-35.5%` at `5 dB` (`40` vs. `62` for `v1` and `160` for `XQRS`)**, lowered `Record 100` discrete `MAE_t` to **`0.27 ms`**, and reduced `5 dB` `RMSSD` bias from **`+2.54%` down to `-0.075%`** via parabolic sub-sample interpolation.

---

## 1. Introduction

Automated QRS detection is foundational to cardiac rhythm monitoring, event synchronization, and Heart Rate Variability (HRV) analysis. Surface ECG potentials (`0.5–3 mV`) are corrupted by respiratory baseline drift (`<0.5 Hz`), `50/60 Hz` power-line interference, electromyographic (`ma`) noise, and mechanical electrode-motion (`em`) artifacts (`1–15 Hz`) that overlap the QRS spectrum directly.

Two methodological blind spots weaken conventional detector benchmarks:
1. **Reporting only `±150 ms` F1-scores**, concealing millisecond-scale R-peak jitter that propagates into `RMSSD` (`RR'_i = RR_i + e_{i+1} - e_i`).
2. **Modifying detector parameters ad hoc without first freezing a baseline and diagnosing physical failure mechanisms.**

This study follows a disciplined three-stage engineering progression: baseline conditioning (`v0.1.0`), failure-mode characterization (`v0.2.0`), and targeted algorithmic improvement (`v0.3.0`).

---

## 2. Materials and Methods

### 2.1 Datasets, Header Integrity & Digital Software Verification
Five 30-minute Channel 0 (`MLII`) MIT-BIH recordings (`100`, `101`, `103`, `105`, `106`; `fs = 360 Hz`, `Δt = 2.78 ms`, `200 counts/mV`, 11-bit ADC) were audited automatically for physical scaling and finiteness, excluding `119` non-beat quality/artifact markers (`~`, `&#124;`) in `Record 105` and `71` rhythm markers (`+`) in `Record 106` (`10,789` true reference beats). Arithmetic scaling and beat-rate routines were verified functionally via a `1.0 mV` / `60.0 bpm` digital pulse train (`0.0%` nominal error; `+5.0%` under `1.05x` gain perturbation).

### 2.2 Baseline (`Pan–Tompkins v1`) vs. Enhanced (`Pan–Tompkins v2`) Pipelines
- **Frozen Baseline (`Pan–Tompkins v1`):** 4th-order zero-phase Butterworth bandpass (`5–20 Hz`, `sosfiltfilt`), centered derivative `d[n]`, squaring, `120 ms` moving-window integration, dual adaptive thresholds (`SPKI`/`NPKI`, `1.66 * mean_RR` searchback, `250 ms` refractory period), and `±90 ms` filtered + `±15 ms` raw-signal peak relocation.
- **Evidence-Driven Upgrade (`Pan–Tompkins v2`):**
  1. *Conditioned Apex Relocation & Parabolic Sub-Sample Vertex:* Removes the noisy `±15 ms` raw-signal snap, centering integer peaks directly on `|y_5_20[n]|` and fitting a 3-point parabola (`r_peaks_subsample`) for continuous HRV analysis.
  2. *Core-Excluded Same-Band `QFSNR` Gate:* Divides peak `5–20 Hz` amplitude (`±60 ms`) by flank RMS (`[-240..-85] U [+85..+240] ms`), rejecting short-RR (`<0.52 * median_RR`) clutter spikes (`QFSNR < 2.8`) and `200–360 ms` motion doublets (`QFSNR < 3.4` with low relative energy) without penalizing wide `130 ms` PVCs.
  3. *Compensatory-Pause-Aware `2.5–15 Hz` Ectopic Rescue:* Estimates clean `sinus_RR` (`0.75–1.35 * median_RR`), separates post-PVC compensatory pauses (`prev_RR < 0.72 * sinus_RR`) from genuine `2-RR` missed-beat gaps (`> 1.56 * sinus_RR`), and rescues low-slope PVCs using a `2.5–15 Hz` zero-phase envelope (`QFSNR >= 3.2`).

---

## 3. Results

1. **Temporal Precision & Downstream `RMSSD` Fidelity (`Record 100`):** While all detectors achieved `100.00%` F1 (`2,265/2,265 TP`, `0 FP`) across synthetic white-noise levels down to `5 dB`, `PT v1` `MAE_t` degraded from `1.50` to `1.96 ms` and `RMSSD` (`n = 2,196` `N–N` pairs) inflated by `+2.54%` (`28.51` vs. `27.80 ms`). `PT v2` cut discrete baseline `MAE_t` to **`0.27 ms`** (`0.47 ms` at `5 dB`, outperforming `WFDB XQRS` at `0.50` and `0.56 ms`) and reduced `5 dB` `RMSSD` bias to **`+0.22%` (discrete)** and **`-0.075%` (parabolic sub-sample, `RR RMSE = 1.35 ms`)**.
2. **Motion & Muscle Artifact Rejection (`Record 105` & `NSTDB`):** `PT v2`'s `QFSNR` gate cut `Record 105` FPs from `33 → 20` at baseline (`F1 = 99.44%`), `34 → 19` at `5 dB` (`F1 = 99.44%`), and `43 → 25` under `+6 dB` physical `nstdb/em` (beating `XQRS`'s `34 FP`), while cutting `-6 dB` physical muscle (`nstdb/ma`) FPs on `Record 100` from `127 → 72` (`F1 = 98.35%` vs. `96.21%` for `XQRS`).
3. **Multiform Ventricular Ectopy (`Record 106`):** Where `WFDB XQRS` missed `88` multiform PVCs at baseline (`Se = 95.65%`, `F1 = 97.77%`) and `93` at `5 dB`, `PT v2` missed only **`1` beat at baseline (`Se = 99.95%`, `F1 = 99.98%`)** and **`3` beats at `5 dB` (`F1 = 99.88%`)**, while raising strict `±8.33 ms` baseline F1 from `95.79%` (`v1`) to **`97.40%`**.
4. **Pooled 5-Record Totals (`10,789` Beats):** Total errors (`FP + FN`) fell from `49` (`PT v1`) and `132` (`XQRS`) to **`36` (`PT v2`)** at baseline, and from `62` (`PT v1`) and `160` (`XQRS`) to **`40` (`PT v2`)** at `5 dB` SNR.

---

## 4. Discussion & Limitations

Freezing a baseline (`v0.2.0`) to isolate exact failure mechanisms before modifying detection logic (`v0.3.0`) resolved the classical sensitivity–specificity trade-off observed between `PT v1` and `WFDB XQRS` on this 5-record benchmark. Remaining limitations include: (1) purposive 5-record pilot scope requiring subject-clustered bootstrap validation across all 48 MIT-BIH records (47 subjects); (2) residual false positives under extreme `-6 dB` electrode motion (`164 FP` on `Record 100`), motivating multi-lead/impedance SQI fusion; and (3) offline zero-phase `sosfiltfilt` conditioning (`5–20 Hz`), which requires causal group-delay redesign for real-time embedded monitors and is not intended for clinical ST/morphology diagnosis.
