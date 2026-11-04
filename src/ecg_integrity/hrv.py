"""Time-domain Heart Rate Variability (HRV) jitter propagation analysis.

Investigates how QRS temporal localization errors (e_i = t_det,i - t_ref,i)
propagate into consecutive normal-to-normal (NN) intervals:
    RR'_i = RR_i + (e_{i+1} - e_i)
and quantify the resulting bias in RMSSD and SDNN under controlled noise.

Methodological safeguard (Section 7):
Conventional time-domain HRV requires clean sinus normal-to-normal ('N' to 'N')
intervals free of ectopic beats ('A', 'V', 'Q', etc.) and missed/false detections.
This module isolates verified consecutive 'N'-'N' intervals that are matched as
True Positives so pure localization jitter is separated from ectopy/beat-loss.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ecg_integrity.evaluation import BeatEvaluationResult


def evaluate_nn_hrv_jitter_impact(
    r_peaks_ref: np.ndarray,
    beat_symbols: list[str],
    eval_result: BeatEvaluationResult,
    fs: float,
    min_rr_ms: float = 300.0,
    max_rr_ms: float = 2000.0,
) -> dict[str, Any]:
    """Compare reference vs. detected NN intervals and RMSSD on verified 'N'-'N' pairs."""
    if not np.isfinite(fs) or fs <= 0.0:
        raise ValueError("fs must be finite and positive.")

    ref_samples = np.asarray(r_peaks_ref, dtype=np.int64)
    if ref_samples.size != len(beat_symbols):
        raise ValueError("r_peaks_ref and beat_symbols must have identical length.")

    # Map original reference beat index -> matched detected sample index
    ref_idx_to_det_sample: dict[int, int] = {
        int(orig_idx): int(det_samp)
        for orig_idx, det_samp in zip(eval_result.matched_ref_indices, eval_result.matched_det_samples)
    }

    nn_ref_ms_list: list[float] = []
    nn_det_ms_list: list[float] = []

    for i in range(ref_samples.size - 1):
        # Strictly require consecutive normal sinus beats ('N' -> 'N')
        if beat_symbols[i] != "N" or beat_symbols[i + 1] != "N":
            continue
        # Require both beats to be matched True Positives to isolate timing-jitter propagation
        if i not in ref_idx_to_det_sample or (i + 1) not in ref_idx_to_det_sample:
            continue

        rr_ref_ms = ((float(ref_samples[i + 1]) - float(ref_samples[i])) / fs) * 1000.0
        if not (min_rr_ms <= rr_ref_ms <= max_rr_ms):
            continue

        rr_det_ms = ((float(ref_idx_to_det_sample[i + 1]) - float(ref_idx_to_det_sample[i])) / fs) * 1000.0
        nn_ref_ms_list.append(rr_ref_ms)
        nn_det_ms_list.append(rr_det_ms)

    nn_ref = np.asarray(nn_ref_ms_list, dtype=np.float64)
    nn_det = np.asarray(nn_det_ms_list, dtype=np.float64)

    if nn_ref.size < 2:
        return {
            "n_nn_intervals": int(nn_ref.size),
            "rmssd_ref_ms": np.nan,
            "rmssd_det_ms": np.nan,
            "rmssd_diff_ms": np.nan,
            "rmssd_error_percent": np.nan,
            "sdnn_ref_ms": np.nan,
            "sdnn_det_ms": np.nan,
            "rr_mae_ms": np.nan,
            "rr_rmse_ms": np.nan,
        }

    rmssd_ref = float(np.sqrt(np.mean(np.diff(nn_ref) ** 2)))
    rmssd_det = float(np.sqrt(np.mean(np.diff(nn_det) ** 2)))
    rmssd_diff = rmssd_det - rmssd_ref
    rmssd_err_pct = ((rmssd_det - rmssd_ref) / rmssd_ref) * 100.0 if rmssd_ref > 0 else np.nan

    sdnn_ref = float(np.std(nn_ref, ddof=1))
    sdnn_det = float(np.std(nn_det, ddof=1))

    rr_errors = nn_det - nn_ref
    rr_mae = float(np.mean(np.abs(rr_errors)))
    rr_rmse = float(np.sqrt(np.mean(rr_errors**2)))

    return {
        "n_nn_intervals": int(nn_ref.size),
        "rmssd_ref_ms": round(rmssd_ref, 4),
        "rmssd_det_ms": round(rmssd_det, 4),
        "rmssd_diff_ms": round(rmssd_diff, 4),
        "rmssd_error_percent": round(float(rmssd_err_pct), 4) if np.isfinite(rmssd_err_pct) else np.nan,
        "sdnn_ref_ms": round(sdnn_ref, 4),
        "sdnn_det_ms": round(sdnn_det, 4),
        "rr_mae_ms": round(rr_mae, 4),
        "rr_rmse_ms": round(rr_rmse, 4),
    }
