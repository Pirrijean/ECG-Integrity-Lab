"""Educational Pan–Tompkins-inspired QRS complex detection algorithm.

Pipeline stages (Step 12):
1. Bandpass filtering (5–20 Hz, 4th-order zero-phase Butterworth SOS)
2. Differentiation: d[n] = dx[n]/dn
3. Squaring: e[n] = d[n]^2
4. Moving-window integration: M[n] = (1/N) * sum_{k=0}^{N-1} e[n - k] (120 ms window)
5. Adaptive Pan–Tompkins dual-threshold peak detection with 250 ms refractory period
6. Local R-peak relocation (±90 ms window) on the filtered waveform with 250 ms minimum separation

Note on temporal resolution:
``detect_qrs_pan_tompkins`` returns integer sample indices ``r_peaks`` whose individual
discrete temporal resolution is ``dt = 1 / fs`` (2.7778 ms at 360 Hz). An optional
independent helper ``refine_peaks_parabolic`` is provided for experimental sub-sample
vertex interpolation without altering the baseline discrete detector.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
from scipy.signal import find_peaks

from ecg_integrity.signal_processing import bandpass_filter, validate_signal_array


@dataclass(frozen=True)
class QRSDetectorConfig:
    """Configurable parameters for the Pan–Tompkins-inspired QRS detector."""

    lowcut_hz: float = 5.0
    highcut_hz: float = 20.0
    filter_order: int = 4
    integration_window_ms: float = 120.0
    refractory_ms: float = 250.0
    search_window_ms: float = 90.0
    min_separation_ms: float = 250.0
    threshold_alpha: float = 0.25
    searchback_rr_factor: float = 1.66

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class QRSDetectionOutput:
    """Detected QRS peak indices, intermediate signal stages, and adaptive threshold trace."""

    r_peaks: np.ndarray
    filtered: np.ndarray
    differentiated: np.ndarray
    squared: np.ndarray
    integrated: np.ndarray
    threshold_trace: np.ndarray
    config: QRSDetectorConfig


def detect_qrs_pan_tompkins(
    signal: np.ndarray,
    fs: float,
    config: QRSDetectorConfig | None = None,
) -> QRSDetectionOutput:
    """Run the complete Pan–Tompkins-inspired QRS detection pipeline."""
    cfg = config or QRSDetectorConfig()
    x = validate_signal_array(signal, fs, max_abs_mv=100.0)

    # 1. Zero-phase Butterworth bandpass (5–20 Hz, order 4)
    filtered = bandpass_filter(
        x,
        fs=fs,
        lowcut_hz=cfg.lowcut_hz,
        highcut_hz=cfg.highcut_hz,
        order=cfg.filter_order,
    )

    # 2. Centered numerical derivative d[n]
    differentiated = np.gradient(filtered)

    # 3. Pointwise squaring e[n] = d[n]^2
    squared = differentiated**2

    # 4. Moving-window integration M[n] over 120 ms window (centered to preserve zero group delay)
    win_samples = max(1, int(round((cfg.integration_window_ms / 1000.0) * fs)))
    kernel = np.ones(win_samples, dtype=np.float64) / float(win_samples)
    integrated = np.convolve(squared, kernel, mode="same")

    # 5. Candidate peak extraction with physiological refractory distance (250 ms)
    refractory_samples = max(1, int(round((cfg.refractory_ms / 1000.0) * fs)))
    candidate_peaks, _ = find_peaks(integrated, distance=refractory_samples)

    if candidate_peaks.size == 0:
        return QRSDetectionOutput(
            r_peaks=np.empty(0, dtype=np.int64),
            filtered=filtered,
            differentiated=differentiated,
            squared=squared,
            integrated=integrated,
            threshold_trace=np.zeros_like(integrated),
            config=cfg,
        )

    # Robust initialization of signal peak (SPKI) and noise peak (NPKI) using 2-second windows
    two_sec = max(1, int(round(2.0 * fs)))
    n_init_windows = min(4, max(1, x.size // two_sec))
    window_maxima = [
        float(np.max(integrated[i * two_sec : min(x.size, (i + 1) * two_sec)]))
        for i in range(n_init_windows)
    ]
    spki = float(np.median(window_maxima))
    npki = float(np.percentile(integrated[: min(x.size, 4 * two_sec)], 50))
    threshold_i1 = npki + cfg.threshold_alpha * (spki - npki)
    threshold_i2 = 0.5 * threshold_i1

    threshold_trace = np.full(x.size, threshold_i1, dtype=np.float64)
    accepted_integrated: list[int] = []
    rr_history: list[int] = []
    rejected_since_last: list[tuple[int, float]] = []
    prev_trace_idx = 0

    for idx in candidate_peaks:
        idx_int = int(idx)
        threshold_trace[prev_trace_idx : idx_int + 1] = threshold_i1
        prev_trace_idx = idx_int + 1
        peak_val = float(integrated[idx_int])

        # Searchback check if long gap exceeds 1.66 * running mean RR
        if accepted_integrated and rr_history and rejected_since_last:
            mean_rr = float(np.mean(rr_history[-8:]))
            if (idx_int - accepted_integrated[-1]) > cfg.searchback_rr_factor * mean_rr:
                best_sb_idx, best_sb_val = max(rejected_since_last, key=lambda item: item[1])
                if best_sb_val >= threshold_i2 and (best_sb_idx - accepted_integrated[-1]) >= refractory_samples:
                    rr_history.append(best_sb_idx - accepted_integrated[-1])
                    accepted_integrated.append(best_sb_idx)
                    spki = 0.25 * min(best_sb_val, 3.0 * spki) + 0.75 * spki
                    rejected_since_last.clear()
                    threshold_i1 = npki + cfg.threshold_alpha * (spki - npki)
                    threshold_i2 = 0.5 * threshold_i1

        if peak_val >= threshold_i1:
            if not accepted_integrated or (idx_int - accepted_integrated[-1]) >= refractory_samples:
                if accepted_integrated:
                    rr_history.append(int(idx_int - accepted_integrated[-1]))
                accepted_integrated.append(idx_int)
                # Clamp extreme transient spikes to prevent post-artifact blindness
                clamped_peak = min(peak_val, 3.0 * spki) if spki > 0 else peak_val
                spki = 0.125 * clamped_peak + 0.875 * spki
                rejected_since_last.clear()
        else:
            npki = 0.125 * peak_val + 0.875 * npki
            rejected_since_last.append((idx_int, peak_val))

        threshold_i1 = npki + cfg.threshold_alpha * (spki - npki)
        threshold_i2 = 0.5 * threshold_i1

    if prev_trace_idx < x.size:
        threshold_trace[prev_trace_idx:] = threshold_i1

    # 6. Relocate each integrated detection onto the local QRS extremum (±90 ms)
    search_samples = max(1, int(round((cfg.search_window_ms / 1000.0) * fs)))
    min_sep_samples = max(1, int(round((cfg.min_separation_ms / 1000.0) * fs)))

    prelim_pos = [
        float(np.max(filtered[max(0, p - search_samples) : min(x.size, p + search_samples + 1)]))
        for p in accepted_integrated
    ]
    prelim_neg = [
        float(np.abs(np.min(filtered[max(0, p - search_samples) : min(x.size, p + search_samples + 1)])))
        for p in accepted_integrated
    ]
    use_positive_polarity = (np.median(prelim_pos) >= 0.85 * np.median(prelim_neg)) if accepted_integrated else True

    relocated: list[int] = []
    for p in accepted_integrated:
        left = max(0, p - search_samples)
        right = min(x.size, p + search_samples + 1)
        seg_filt = filtered[left:right]
        if use_positive_polarity and float(np.max(seg_filt)) > 0.4 * float(np.max(np.abs(seg_filt))):
            local_offset = int(np.argmax(seg_filt))
        else:
            local_offset = int(np.argmax(np.abs(seg_filt)))
        cand = left + local_offset

        # Fine local snap (±15 ms) onto the band-preserved raw deviation from local median
        fine_radius = max(1, int(round(0.015 * fs)))
        f_left = max(0, cand - fine_radius)
        f_right = min(x.size, cand + fine_radius + 1)
        local_baseline = float(np.median(x[max(0, cand - search_samples) : min(x.size, cand + search_samples + 1)]))
        seg_raw_dev = x[f_left:f_right] - local_baseline
        if use_positive_polarity and float(np.max(seg_raw_dev)) > 0.4 * float(np.max(np.abs(seg_raw_dev))):
            cand = f_left + int(np.argmax(seg_raw_dev))
        else:
            cand = f_left + int(np.argmax(np.abs(seg_raw_dev)))

        if relocated and (cand - relocated[-1]) < min_sep_samples:
            if abs(filtered[cand]) > abs(filtered[relocated[-1]]):
                relocated[-1] = cand
        else:
            relocated.append(cand)

    return QRSDetectionOutput(
        r_peaks=np.asarray(relocated, dtype=np.int64),
        filtered=filtered,
        differentiated=differentiated,
        squared=squared,
        integrated=integrated,
        threshold_trace=threshold_trace,
        config=cfg,
    )


def detect_qrs(
    signal: np.ndarray,
    fs: float,
    config: QRSDetectorConfig | None = None,
) -> np.ndarray:
    """Convenience wrapper returning only the detected R-peak integer sample indices."""
    return detect_qrs_pan_tompkins(signal, fs=fs, config=config).r_peaks


def refine_peaks_parabolic(
    waveform: np.ndarray,
    peak_indices: np.ndarray,
) -> np.ndarray:
    """Optional sub-sample 3-point parabolic vertex interpolation around discrete peaks.

    Provided for independent methodological comparison when sub-sample fiducial
    refinement is studied; NOT enabled in the baseline discrete detector.
    """
    arr = np.asarray(waveform, dtype=np.float64)
    peaks = np.asarray(peak_indices, dtype=np.int64)
    refined = peaks.astype(np.float64).copy()

    for idx, k in enumerate(peaks):
        if 0 < k < (arr.size - 1):
            y_minus = abs(arr[k - 1])
            y_zero = abs(arr[k])
            y_plus = abs(arr[k + 1])
            denom = y_minus - 2.0 * y_zero + y_plus
            if abs(denom) > 1e-12 and y_zero >= y_minus and y_zero >= y_plus:
                delta = 0.5 * (y_minus - y_plus) / denom
                if abs(delta) <= 0.5:
                    refined[idx] = float(k) + delta
    return refined


@dataclass(frozen=True)
class QRSDetectorV2Config:
    """Configuration for Pan–Tompkins v2 (v0.3.0 evidence-driven enhanced detector).

    Addresses the three failure mechanisms characterized in v0.2.0:
    1. Conditioned 5–20 Hz R-peak relocation + optional parabolic sub-sample vertex
       (eliminating unfiltered high-frequency raw-sample snap bias and AWGN jitter).
    2. Same-band QRS-to-Flank Signal Quality Index (QFSNR, excluding ±85 ms QRS core)
       post-screening gate to reject interpolated electrode-motion / EMG artifacts.
    3. Compensatory-pause-aware ectopic envelope rescue (2.5–15 Hz) + optional
       secondary lead fusion to recover low-slope multiform PVCs.
    """

    lowcut_hz: float = 5.0
    highcut_hz: float = 20.0
    filter_order: int = 4
    ectopic_lowcut_hz: float = 2.5
    ectopic_highcut_hz: float = 15.0
    ectopic_filter_order: int = 3
    integration_window_ms: float = 120.0
    refractory_ms: float = 250.0
    search_window_ms: float = 85.0
    min_separation_ms: float = 250.0
    threshold_alpha: float = 0.25
    searchback_rr_factor: float = 1.66
    enable_sqi_artifact_gate: bool = True
    enable_ectopic_rescue: bool = True

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class QRSDetectionV2Output:
    """Output bundle for Pan–Tompkins v2 including discrete & sub-sample R-peaks and SQI diagnostics."""

    r_peaks: np.ndarray
    r_peaks_subsample: np.ndarray
    filtered: np.ndarray
    ectopic_filtered: np.ndarray
    differentiated: np.ndarray
    squared: np.ndarray
    integrated: np.ndarray
    ectopic_integrated: np.ndarray
    threshold_trace: np.ndarray
    sqi_rejected_count: int
    ectopic_rescued_count: int
    config: QRSDetectorV2Config


def detect_qrs_pan_tompkins_v2(
    signal: np.ndarray,
    fs: float,
    config: QRSDetectorV2Config | None = None,
    secondary_signal: np.ndarray | None = None,
) -> QRSDetectionV2Output:
    """Run the enhanced Pan–Tompkins v2 QRS detector (v0.3.0).

    Preserves ``detect_qrs_pan_tompkins`` (v1) unchanged for historical baseline
    reproducibility while integrating three empirically validated improvements:
    - Stage 5b: QRS-to-Flank SNR (QFSNR) artifact screening for interpolated motion/EMG transients.
    - Stage 5c: Compensatory-pause-aware 2.5–15 Hz ectopic envelope rescue for low-dV/dt PVCs.
    - Stage 6: Zero-phase 5–20 Hz conditioned apex relocation + 3-point parabolic sub-sample output.
    """
    cfg = config or QRSDetectorV2Config()
    x = validate_signal_array(signal, fs, max_abs_mv=100.0)

    # 1. Zero-phase Butterworth bands: QRS (5–20 Hz), wide ectopic (2.5–15 Hz), and <2 Hz drift
    filtered = bandpass_filter(
        x,
        fs=fs,
        lowcut_hz=cfg.lowcut_hz,
        highcut_hz=cfg.highcut_hz,
        order=cfg.filter_order,
    )
    filt_ect = bandpass_filter(
        x,
        fs=fs,
        lowcut_hz=cfg.ectopic_lowcut_hz,
        highcut_hz=cfg.ectopic_highcut_hz,
        order=cfg.ectopic_filter_order,
    )
    nyq_safe = min(80.0, 0.45 * float(fs))
    lf_drift = x - bandpass_filter(x, fs=fs, lowcut_hz=2.0, highcut_hz=nyq_safe, order=3)

    # 2. Derivative & squaring (with optional normalized dual-lead fusion if secondary_signal is given)
    differentiated = np.gradient(filtered)
    squared = differentiated**2
    if secondary_signal is not None:
        x_sec = validate_signal_array(secondary_signal, fs, max_abs_mv=100.0)
        if x_sec.size != x.size:
            raise ValueError(f"secondary_signal length ({x_sec.size}) must equal primary signal length ({x.size}).")
        sec_filt = bandpass_filter(x_sec, fs=fs, lowcut_hz=cfg.lowcut_hz, highcut_hz=cfg.highcut_hz, order=cfg.filter_order)
        sec_sq = np.gradient(sec_filt) ** 2
        s0 = float(np.percentile(squared, 96)) + 1e-12
        s1 = float(np.percentile(sec_sq, 96)) + 1e-12
        squared = (0.65 * (squared / s0) + 0.35 * (sec_sq / s1)) * s0

    # 3. Moving-window integration on primary derivative energy and ectopic amplitude energy
    win_samples = max(1, int(round((cfg.integration_window_ms / 1000.0) * fs)))
    kernel = np.ones(win_samples, dtype=np.float64) / float(win_samples)
    integrated = np.convolve(squared, kernel, mode="same")
    ectopic_integrated = np.convolve(filt_ect**2, kernel, mode="same")

    refractory_samples = max(1, int(round((cfg.refractory_ms / 1000.0) * fs)))
    candidate_peaks, _ = find_peaks(integrated, distance=refractory_samples)
    ectopic_peaks, _ = find_peaks(ectopic_integrated, distance=refractory_samples)

    if candidate_peaks.size == 0:
        empty_i = np.empty(0, dtype=np.int64)
        return QRSDetectionV2Output(
            r_peaks=empty_i,
            r_peaks_subsample=np.empty(0, dtype=np.float64),
            filtered=filtered,
            ectopic_filtered=filt_ect,
            differentiated=differentiated,
            squared=squared,
            integrated=integrated,
            ectopic_integrated=ectopic_integrated,
            threshold_trace=np.zeros_like(integrated),
            sqi_rejected_count=0,
            ectopic_rescued_count=0,
            config=cfg,
        )

    # 4. Adaptive dual-threshold Pan–Tompkins primary pass
    two_sec = max(1, int(round(2.0 * fs)))
    n_init_windows = min(4, max(1, x.size // two_sec))
    spki = float(
        np.median([float(np.max(integrated[i * two_sec : min(x.size, (i + 1) * two_sec)])) for i in range(n_init_windows)])
    )
    npki = float(np.percentile(integrated[: min(x.size, 4 * two_sec)], 50))
    ect_spki = float(
        np.median(
            [float(np.max(ectopic_integrated[i * two_sec : min(x.size, (i + 1) * two_sec)])) for i in range(n_init_windows)]
        )
    )

    threshold_i1 = npki + cfg.threshold_alpha * (spki - npki)
    threshold_i2 = 0.5 * threshold_i1
    threshold_trace = np.full(x.size, threshold_i1, dtype=np.float64)

    accepted: list[int] = []
    rr_history: list[int] = []
    rejected_since_last: list[tuple[int, float]] = []
    prev_trace_idx = 0

    for idx in candidate_peaks:
        idx_int = int(idx)
        threshold_trace[prev_trace_idx : idx_int + 1] = threshold_i1
        prev_trace_idx = idx_int + 1
        peak_val = float(integrated[idx_int])

        if accepted and rr_history and rejected_since_last:
            mean_rr = float(np.mean(rr_history[-8:]))
            if (idx_int - accepted[-1]) > cfg.searchback_rr_factor * mean_rr:
                best_sb_idx, best_sb_val = max(rejected_since_last, key=lambda item: item[1])
                if best_sb_val >= threshold_i2 and (best_sb_idx - accepted[-1]) >= refractory_samples:
                    rr_history.append(best_sb_idx - accepted[-1])
                    accepted.append(best_sb_idx)
                    spki = 0.25 * min(best_sb_val, 3.0 * spki) + 0.75 * spki
                    rejected_since_last.clear()
                    threshold_i1 = npki + cfg.threshold_alpha * (spki - npki)
                    threshold_i2 = 0.5 * threshold_i1

        if peak_val >= threshold_i1:
            if not accepted or (idx_int - accepted[-1]) >= refractory_samples:
                if accepted:
                    rr_history.append(int(idx_int - accepted[-1]))
                accepted.append(idx_int)
                clamped_peak = min(peak_val, 3.0 * spki) if spki > 0 else peak_val
                spki = 0.125 * clamped_peak + 0.875 * spki
                rejected_since_last.clear()
        else:
            npki = 0.125 * peak_val + 0.875 * npki
            rejected_since_last.append((idx_int, peak_val))

        threshold_i1 = npki + cfg.threshold_alpha * (spki - npki)
        threshold_i2 = 0.5 * threshold_i1

    if prev_trace_idx < x.size:
        threshold_trace[prev_trace_idx:] = threshold_i1

    # Local same-band QRS-to-Flank SNR helper (excludes ±85 ms QRS core so wide PVCs are never self-penalized)
    w60 = max(1, int(round(0.060 * fs)))
    w85 = max(1, int(round(0.085 * fs)))
    w240 = max(1, int(round(0.240 * fs)))

    def compute_qfsnr(c: int) -> tuple[float, float, float, float]:
        l60, r60 = max(0, c - w60), min(x.size, c + w60 + 1)
        pk5 = float(np.max(np.abs(filtered[l60:r60])))
        fl5 = np.concatenate([filtered[max(0, c - w240) : max(0, c - w85)], filtered[min(x.size, c + w85) : min(x.size, c + w240)]])
        fl5_rms = float(np.sqrt(np.mean(fl5**2))) + 1e-9 if fl5.size > 0 else 1e-9

        pke = float(np.max(np.abs(filt_ect[l60:r60])))
        fle = np.concatenate([filt_ect[max(0, c - w240) : max(0, c - w85)], filt_ect[min(x.size, c + w85) : min(x.size, c + w240)]])
        fle_rms = float(np.sqrt(np.mean(fle**2))) + 1e-9 if fle.size > 0 else 1e-9

        pre_lf = lf_drift[max(0, c - w240) : max(1, c - w85)]
        post_lf = lf_drift[min(x.size - 1, c + w85) : min(x.size, c + w240)]
        dc_jump = abs(float(np.median(pre_lf)) - float(np.median(post_lf))) if (pre_lf.size > 0 and post_lf.size > 0) else 0.0
        raw_ptp = float(np.ptp(x[max(0, c - w240) : min(x.size, c + w240)]))
        return pk5 / fl5_rms, pke / fle_rms, dc_jump, raw_ptp

    # 5b. Signal Quality Index (QFSNR) Artifact Screening Gate
    sqi_rejected_count = 0
    if cfg.enable_sqi_artifact_gate and len(accepted) >= 6:
        rr_med = float(np.median(np.diff(accepted)))
        cleaned: list[int] = []
        for i, c in enumerate(accepted):
            if cleaned:
                prev_c = cleaned[-1]
                dt_rr = (c - prev_c) / rr_med
                if dt_rr < 0.72 and (i + 1) < len(accepted):
                    span = (accepted[i + 1] - prev_c) / rr_med
                    qf5_p, _, _, _ = compute_qfsnr(prev_c)
                    qf5_c, qfe_c, jmp_c, ptp_c = compute_qfsnr(c)
                    # Case 1: Interpolated motion/EMG transient splitting one normal sinus cycle (~1.0 * RR)
                    if 0.78 <= span <= 1.26:
                        if (
                            qf5_c < 2.65
                            or (qf5_c < 3.45 and (jmp_c > 0.75 or ptp_c > 3.40))
                            or (qf5_c < 0.45 * qf5_p and qfe_c < 2.65)
                        ):
                            sqi_rejected_count += 1
                            continue
                    # Case 2: Ultra-close (<0.48 * RR) artifact doublet preceding/following a clean QRS
                    if dt_rr < 0.48 and qf5_p < 2.35 and qf5_c >= 3.40:
                        cleaned[-1] = c
                        sqi_rejected_count += 1
                        continue
                    if dt_rr < 0.48 and qf5_c < 2.35 and qf5_p >= 3.40:
                        sqi_rejected_count += 1
                        continue
            cleaned.append(c)
        accepted = cleaned

    # 5c. Compensatory-Pause-Aware Ectopic Envelope Rescue (2.5–15 Hz) for low-slope multiform PVCs
    ectopic_rescued_count = 0
    if cfg.enable_ectopic_rescue and len(accepted) >= 6:
        all_rr = np.diff(accepted)
        m0 = float(np.median(all_rr))
        sinus_rrs = all_rr[(all_rr >= 0.75 * m0) & (all_rr <= 1.35 * m0)]
        sinus_rr = float(np.median(sinus_rrs)) if sinus_rrs.size >= 3 else m0
        margin = max(refractory_samples, int(round(0.28 * fs)))

        rescued: list[int] = [accepted[0]]
        for k_i in range(1, len(accepted)):
            prev_c = rescued[-1]
            cur_c = accepted[k_i]
            gap = cur_c - prev_c
            prev_rr = (prev_c - rescued[-2]) if len(rescued) >= 2 else sinus_rr
            # Distinguish physiological compensatory pause after a premature beat from a 2-RR missed-beat gap
            is_comp_pause = (prev_rr < 0.72 * sinus_rr) and (gap < 1.62 * sinus_rr)
            if gap > 1.56 * sinus_rr and not is_comp_pause:
                lo = prev_c + margin
                hi = cur_c - margin
                if hi > lo:
                    ec_in = ectopic_peaks[(ectopic_peaks >= lo) & (ectopic_peaks <= hi)]
                    best_c: int | None = None
                    best_v = -1.0
                    for ec in ec_in:
                        eci = int(ec)
                        _, qfe, jmp, _ = compute_qfsnr(eci)
                        ev = float(ectopic_integrated[eci])
                        if qfe >= 2.20 and jmp < 1.50 and ev >= 0.12 * ect_spki and ev > best_v:
                            best_c = eci
                            best_v = ev
                    if best_c is not None:
                        rescued.append(best_c)
                        ectopic_rescued_count += 1
            rescued.append(cur_c)
        accepted = rescued

    # 6. Zero-phase 5–20 Hz conditioned R-peak relocation (eliminates unfiltered HF raw-snap bias & jitter)
    search_samples = max(1, int(round((cfg.search_window_ms / 1000.0) * fs)))
    min_sep_samples = max(1, int(round((cfg.min_separation_ms / 1000.0) * fs)))

    prelim_pos = [
        float(np.max(filtered[max(0, p - search_samples) : min(x.size, p + search_samples + 1)]))
        for p in accepted
    ]
    prelim_neg = [
        float(np.abs(np.min(filtered[max(0, p - search_samples) : min(x.size, p + search_samples + 1)])))
        for p in accepted
    ]
    use_positive_polarity = (np.median(prelim_pos) >= 0.85 * np.median(prelim_neg)) if accepted else True

    relocated: list[int] = []
    for p in accepted:
        left = max(0, p - search_samples)
        right = min(x.size, p + search_samples + 1)
        seg_filt = filtered[left:right]
        if use_positive_polarity and float(np.max(seg_filt)) >= 0.45 * float(np.max(np.abs(seg_filt))):
            cand = left + int(np.argmax(seg_filt))
        else:
            cand = left + int(np.argmax(np.abs(seg_filt)))

        if relocated and (cand - relocated[-1]) < min_sep_samples:
            if abs(filtered[cand]) > abs(filtered[relocated[-1]]):
                relocated[-1] = cand
        else:
            relocated.append(cand)

    r_peaks_arr = np.asarray(relocated, dtype=np.int64)
    r_peaks_sub = refine_peaks_parabolic(filtered, r_peaks_arr)

    return QRSDetectionV2Output(
        r_peaks=r_peaks_arr,
        r_peaks_subsample=r_peaks_sub,
        filtered=filtered,
        ectopic_filtered=filt_ect,
        differentiated=differentiated,
        squared=squared,
        integrated=integrated,
        ectopic_integrated=ectopic_integrated,
        threshold_trace=threshold_trace,
        sqi_rejected_count=sqi_rejected_count,
        ectopic_rescued_count=ectopic_rescued_count,
        config=cfg,
    )
