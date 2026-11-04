# Beyond F1-Score: Evaluating Temporal Accuracy and Noise Sensitivity of a Pan–Tompkins-Inspired QRS Detection Pipeline

*An exploratory evaluation using the MIT-BIH Arrhythmia and Noise Stress Test Databases*

**Author:** Pirrijean  
**Repository:** [https://github.com/Pirrijean/ECG-Integrity-Lab](https://github.com/Pirrijean/ECG-Integrity-Lab) (`v0.2.0`)  

---

## Abstract

Standard evaluations of electrocardiogram (ECG) QRS detectors routinely report aggregate beat-classification metrics (Sensitivity, Positive Predictive Value, and F1-score) under a single permissive matching window (typically $\pm 150\text{ ms}$). Using **ECG-Integrity-Lab**, an open Python verification framework, we evaluate a zero-phase Butterworth ($5\text{–}20\text{ Hz}$) and Pan–Tompkins-inspired QRS detection pipeline across five full-duration (30-minute) recordings from the MIT-BIH Arrhythmia Database (`10,789` expert-reviewed reference beats), physical ambulatory noise from the MIT-BIH Noise Stress Test Database (`nstdb`: `bw`, `ma`, `em`), six temporal matching tolerances ($8.33\text{–}150\text{ ms}$), and an independent benchmark against WFDB's `XQRS` detector. Pooled baseline counts yielded `99.89%` Sensitivity (`10,777 TP`, `12 FN`) and `99.66%` PPV (`37 FP`), yet errors were strongly concentrated by signal quality and beat morphology (`33/37` FPs in motion-corrupted `Record 105`). Crucially, we show that **classification metrics can remain at 100.00% while temporal localization accuracy degrades substantially**: on `Record 100`, decreasing white-noise SNR from baseline to $5\text{ dB}$ produced zero false positives or missed beats, but increased mean absolute R-peak timing error ($\text{MAE}_t$) by $+30.7\%$ ($1.50 \rightarrow 1.96\text{ ms}$ across a $2.78\text{ ms}$ discrete sampling interval), timing error standard deviation by $+56.0\%$ ($1.68 \rightarrow 2.62\text{ ms}$), and normal-to-normal `RMSSD` heart-rate variability estimates by $+2.54\%$. Tightening the matching tolerance from $\pm 150\text{ ms}$ to $\pm 8.33\text{ ms}$ reduced `Record 106` baseline F1 from `99.95%` to `95.79%` due to fiducial dispersion on multiform premature ventricular contractions (PVCs), while physical electrode motion (`nstdb/em` at $-6\text{ dB}$) induced `171 FP` (`PPV = 92.97%`), demonstrating why synthetic sinusoidal perturbations alone cannot substantiate claims of ambulatory noise robustness.

---

## 1. Introduction

Automated QRS complex detection underpins heart-rate monitoring, physiological event synchronization, and heart-rate variability (HRV) analysis. Surface ECG potentials ($0.5\text{–}3\text{ mV}$) are subject to respiratory baseline drift ($<0.5\text{ Hz}$), $50/60\text{ Hz}$ power-line coupling, electromyographic (EMG) interference, and mechanical electrode-motion transients ($1\text{–}15\text{ Hz}$) whose spectral energy overlaps the QRS complex directly.

Two methodological blind spots frequently weaken engineering evaluations of QRS algorithms:
1. **Reliance on a single wide tolerance ($\pm 150\text{ ms}$) and count-based F1-score**, which conceals temporal R-peak jitter relevant to HRV and waveform alignment.
2. **Equating synthetic sinusoidal/Gaussian perturbations with real ambulatory acquisition artifacts.**

This study characterizes both beat-detection performance and temporal localization precision under controlled synthetic perturbations, physical `nstdb` noise recordings, morphology-stratified error diagnostics, and multi-tolerance benchmarking.

---

## 2. Materials and Methods

### 2.1 Datasets and Pre-Processing Integrity Audit
Five 30-minute Channel 0 (`MLII`) recordings (`100`, `101`, `103`, `105`, `106`) from the MIT-BIH Arrhythmia Database ($f_s = 360\text{ Hz}$, $\Delta t = 2.78\text{ ms}$, `200 counts/mV`, 11-bit ADC) were selected purposively to span clean sinus rhythm, baseline drift, severe tape/motion noise (`105`), and frequent multiform PVCs (`106`, $25.6\%$ `V` beats). Header units (`mV`), calibration gains, and sample finiteness were verified automatically prior to processing. Expert `.atr` annotations were filtered strictly to standard PhysioNet QRS beat symbols, excluding `119` non-beat quality/artifact markers (`~`, `|`) in `Record 105` and `71` rhythm markers (`+`) in `Record 106`.

### 2.2 Functional Software Verification, Conditioning, and Detection
Arithmetic scaling and interval routines were verified functionally using a synthetic $1.0\text{ mV}$ / $60.0\text{ bpm}$ digital pulse train ($0.0\%$ nominal amplitude/rate error; $+5.0\%$ measured amplitude error under a $1.05\times$ scale perturbation). Signal conditioning applied a 4th-order Butterworth bandpass ($5\text{–}20\text{ Hz}$) in second-order sections via zero-phase `sosfiltfilt` (scoped strictly to QRS enhancement, not diagnostic ST/P/T preservation), followed by centered differentiation $d[n]$, squaring $e[n] = d[n]^2$, $120\text{ ms}$ moving-window integration $M[n]$, dual Pan–Tompkins adaptive thresholds (`SPKI`/`NPKI` with $1.66\times\overline{RR}$ searchback and $250\text{ ms}$ refractory period), and $\pm 90\text{ ms}$ local peak relocation.

### 2.3 Perturbation Protocols and Evaluation Metrics
Each recording was evaluated under (i) synthetic baseline wander ($0.20\text{ mV}$, $0.25\text{ Hz}$), $60\text{ Hz}$ interference ($0.05\text{ mV}$), and zero-mean AWGN ($20, 10, 5\text{ dB}$ SNR referenced to $\text{Var}(x)$), and (ii) physical `nstdb` recordings (`bw`, `ma`, `em`) scaled via $\alpha = \sqrt{P_x / (P_n \cdot 10^{\text{SNR}/10})}$ from $24$ down to $-6\text{ dB}$ SNR. Strict one-to-one beat matching was swept across $\pm\{8.33, 16.67, 25, 50, 100, 150\}\text{ ms}$. Signed timing errors followed $e_i = t_{\text{det},i} - t_{\text{ref},i}$. Because the `10,789` beats evaluated across 6 conditions represent `64,734` paired repeated-measures opportunities rather than independent observations, baseline summaries were computed by pooling counts (`compute_micro_aggregated_metrics`).

---

## 3. Results

1. **Pooled Baseline & Tolerance Sensitivity:** Pooled baseline performance at $\pm 150\text{ ms}$ reached `Se = 99.89%`, `PPV = 99.66%`, and `F1 = 99.77%` (`10,777 TP`, `37 FP`, `12 FN`). However, shrinking the tolerance window exposed substantial fiducial dispersion on ectopic beats: on `Record 106`, baseline F1 remained `99.95%` at $\pm 150\text{ ms}$ and $\pm 100\text{ ms}$, `99.80%` at $\pm 50\text{ ms}$, `97.38%` at $\pm 25\text{ ms}$, and dropped to **`95.79%`** (`93.96%` at $5\text{ dB}$ AWGN) at $\pm 8.33\text{ ms}$ ($\pm 3$ samples).
2. **Temporal Jitter & HRV Propagation at Constant 100% F1:** On `Record 100` (30 min, $n=2,265$), F1 stayed at `100.00%` across all synthetic conditions, whereas ensemble $\text{MAE}_t$ rose from `1.50 ms` to `1.96 ms` ($+30.7\%$) across the $2.78\text{ ms}$ discrete sampling grid, $\sigma_t$ increased from `1.68 ms` to `2.62 ms` ($+56.0\%$), and `RMSSD` across `2,196` clean `N`–`N` sinus intervals increased from `27.80 ms` (reference) to `28.51 ms` ($+2.54\%$) at $5\text{ dB}$ SNR.
3. **Synthetic vs. Real `NSTDB` Noise:** Stationary $0.25\text{ Hz}$ and $60\text{ Hz}$ sinusoids and physical `nstdb/bw` down to $-6\text{ dB}$ caused zero added FPs on `Record 100`. Conversely, physical EMG (`nstdb/ma`) and electrode motion (`nstdb/em`) at $-6\text{ dB}$ caused `127 FP` (`PPV = 94.68%`) and `171 FP` (`PPV = 92.97%`, $\sigma_t = 6.38\text{ ms}$), respectively.
4. **Comparison with WFDB `XQRS`:** On `Record 100` and `Record 105`, `XQRS` achieved lower $\text{MAE}_t$ (`0.50 ms` vs. `1.50 ms` on `100`) and fewer FPs under `nstdb/em` (`34` vs. `43` at $+6\text{ dB}$ on `105`). Conversely, on `Record 106` (`518` multiform PVCs), `XQRS` missed `88` beats at baseline (`Se = 95.65%`, `F1 = 97.77%`), whereas the derivative-energy Pan–Tompkins pipeline missed only `2` (`Se = 99.90%`, `F1 = 99.95%`).

---

## 4. Error Analysis (`Records 105 & 106`)

- **Motion-Artifact False Positives (`Record 105`):** Waveform inspection ($\pm 300\text{ ms}$) of all `33` baseline FPs showed that `3` coincided within $3\text{–}36\text{ ms}$ of PhysioNet `|` isolated-artifact annotations, and `19` lay within $5\text{ s}$ of `~` signal-quality transitions (`raw_ptp_mv` up to `5.77 mV`), producing in-band $5\text{–}20\text{ Hz}$ energy bursts up to $27.7\times$ `THRESHOLD_I1`.
- **Biophysical PVC Attenuation (`Record 106` at $5\text{ dB}$):** Normal `N` beats achieved `100.00%` sensitivity (`1503/1503`), whereas `6/518` multiform PVCs were missed (`Se = 98.84%`). Missed PVCs exhibited lower peak-to-peak amplitude ($1.50\text{ mV}$ vs. $2.57\text{ mV}$), $3.2\times$ lower filtered slope ($18.6\text{ mV/s}$ vs. $58.8\text{ mV/s}$), and integrated energy ratios strictly below threshold ($E_{\text{QRS}} / I_1 < 0.60$ for all `6` FNs vs. $2.89$ mean for `512` TPs).

---

## 5. Discussion & Limitations

These findings demonstrate that evaluating QRS detectors solely via $\pm 150\text{ ms}$ F1-scores under synthetic Gaussian/sinusoidal noise masks clinically and technically meaningful failure modes. Key limitations of the present framework include: (1) purposive 5-record pilot scope rather than full 48-record / 47-subject clustered bootstrap inference; (2) single-lead (`MLII`) offline non-causal filtering (`sosfiltfilt`); (3) non-diagnostic $5\text{–}20\text{ Hz}$ passband unsuitable for ST/morphology analysis; and (4) discrete sample-index quantization ($\Delta t = 2.78\text{ ms}$), motivating independent validation of sub-sample interpolation methods.

---

## 6. Conclusion

`ECG-Integrity-Lab` establishes a reproducible methodology for characterizing both beat-detection errors and sub-tolerance temporal localization drift. Reporting tolerance-sweep curves, physical `nstdb` stress tests, morphology-stratified metrics, and NN-interval `RMSSD` sensitivity provides a far more rigorous engineering assessment than single-threshold classification accuracy alone.
