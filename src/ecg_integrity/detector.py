"""Educational Pan–Tompkins-inspired QRS complex detection algorithm.

Pipeline stages (Step 12):
1. Bandpass filtering (5–20 Hz, 4th-order zero-phase Butterworth SOS)
2. Differentiation: d[n] = dx[n]/dn
3. Squaring: e[n] = d[n]^2
4. Moving-window integration: M[n] = (1/N) * sum_{k=0}^{N-1} e[n - k] (120 ms window)
5. Adaptive Pan–Tompkins dual-threshold peak detection with 250 ms refractory period
6. Local R-peak relocation (±90 ms window) on the filtered waveform with 250 ms minimum separation
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
    """Detected QRS peak indices and intermediate signal stages."""

    r_peaks: np.ndarray
    filtered: np.ndarray
    differentiated: np.ndarray
    squared: np.ndarray
    integrated: np.ndarray
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

    accepted_integrated: list[int] = []
    rr_history: list[int] = []
    rejected_since_last: list[tuple[int, float]] = []

    for idx in candidate_peaks:
        peak_val = float(integrated[idx])

        # Searchback check if long gap exceeds 1.66 * running mean RR
        if accepted_integrated and rr_history and rejected_since_last:
            mean_rr = float(np.mean(rr_history[-8:]))
            if (idx - accepted_integrated[-1]) > cfg.searchback_rr_factor * mean_rr:
                best_sb_idx, best_sb_val = max(rejected_since_last, key=lambda item: item[1])
                if best_sb_val >= threshold_i2 and (best_sb_idx - accepted_integrated[-1]) >= refractory_samples:
                    rr_history.append(best_sb_idx - accepted_integrated[-1])
                    accepted_integrated.append(best_sb_idx)
                    spki = 0.25 * min(best_sb_val, 3.0 * spki) + 0.75 * spki
                    rejected_since_last.clear()
                    threshold_i1 = npki + cfg.threshold_alpha * (spki - npki)
                    threshold_i2 = 0.5 * threshold_i1

        if peak_val >= threshold_i1:
            if not accepted_integrated or (idx - accepted_integrated[-1]) >= refractory_samples:
                if accepted_integrated:
                    rr_history.append(int(idx - accepted_integrated[-1]))
                accepted_integrated.append(int(idx))
                # Clamp extreme transient spikes to prevent post-artifact blindness
                clamped_peak = min(peak_val, 3.0 * spki) if spki > 0 else peak_val
                spki = 0.125 * clamped_peak + 0.875 * spki
                rejected_since_last.clear()
        else:
            npki = 0.125 * peak_val + 0.875 * npki
            rejected_since_last.append((int(idx), peak_val))

        threshold_i1 = npki + cfg.threshold_alpha * (spki - npki)
        threshold_i2 = 0.5 * threshold_i1

    # 6. Relocate each integrated detection onto the local QRS extremum (±90 ms)
    search_samples = max(1, int(round((cfg.search_window_ms / 1000.0) * fs)))
    min_sep_samples = max(1, int(round((cfg.min_separation_ms / 1000.0) * fs)))

    # Determine dominant QRS deflection polarity on filtered signal so relocation snaps consistently
    # to the primary R/S wave rather than switching between positive R and negative S lobes.
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
        config=cfg,
    )


def detect_qrs(
    signal: np.ndarray,
    fs: float,
    config: QRSDetectorConfig | None = None,
) -> np.ndarray:
    """Convenience wrapper returning only the detected R-peak sample indices."""
    return detect_qrs_pan_tompkins(signal, fs=fs, config=config).r_peaks
