"""Quantitative QRS detection evaluator with strict one-to-one beat matching.

Implements:
- Strict one-to-one temporal matching within a declared tolerance window.
- Sensitivity (Se), Positive Predictive Value (PPV), F1-score, and Mean Absolute
  Temporal Error (MAE_t in ms).
- Multi-window temporal tolerance sweep (8.33, 16.67, 25, 50, 100, 150 ms).
- Micro-aggregated count summarization across multiple records.

Important methodological conventions:
1. Signed timing error convention:
       e_i = t_det,i - t_ref,i = (sample_det,i - sample_ref,i) / fs
   A positive value indicates that the detected peak occurs AFTER the expert
   reference annotation; a negative value indicates it precedes the reference.
2. Discrete sampling interval vs. ensemble MAE_t:
   At fs = 360 Hz, individual integer-sample detections have a discrete resolution
   of Δt = 1 / 360 s ≈ 2.7778 ms. Values of MAE_t < 2.78 ms represent the arithmetic
   mean of discrete errors (including 0-sample offsets) across many beats, NOT
   sub-sample resolution of single detections.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

# Standard temporal matching tolerances (ms) corresponding to 3, 6, 9, 18, 36, and 54 samples at 360 Hz
DEFAULT_TOLERANCES_MS: tuple[float, ...] = (8.33, 16.67, 25.0, 50.0, 100.0, 150.0)


@dataclass(frozen=True)
class BeatEvaluationResult:
    """Complete quantitative evaluation metrics for a QRS detection run."""

    n_ref: int
    n_det: int
    tp: int
    fp: int
    fn: int
    sensitivity: float
    ppv: float
    f1_score: float
    mae_ms: float
    mean_signed_error_ms: float
    std_error_ms: float
    median_abs_error_ms: float
    max_abs_error_ms: float
    tolerance_ms: float
    sampling_interval_ms: float
    matched_ref_samples: np.ndarray
    matched_det_samples: np.ndarray
    matched_ref_indices: np.ndarray
    timing_errors_ms: np.ndarray
    unmatched_ref_samples: np.ndarray
    unmatched_ref_indices: np.ndarray
    unmatched_det_samples: np.ndarray

    def to_metrics_dict(self) -> dict[str, Any]:
        """Return scalar summary metrics suitable for CSV/JSON serialization."""
        return {
            "n_ref": self.n_ref,
            "n_det": self.n_det,
            "tp": self.tp,
            "fp": self.fp,
            "fn": self.fn,
            "sensitivity": round(self.sensitivity, 6),
            "ppv": round(self.ppv, 6),
            "f1_score": round(self.f1_score, 6),
            "mae_ms": round(self.mae_ms, 4) if np.isfinite(self.mae_ms) else np.nan,
            "mean_signed_error_ms": (
                round(self.mean_signed_error_ms, 4) if np.isfinite(self.mean_signed_error_ms) else np.nan
            ),
            "std_error_ms": round(self.std_error_ms, 4) if np.isfinite(self.std_error_ms) else np.nan,
            "median_abs_error_ms": (
                round(self.median_abs_error_ms, 4) if np.isfinite(self.median_abs_error_ms) else np.nan
            ),
            "max_abs_error_ms": (
                round(self.max_abs_error_ms, 4) if np.isfinite(self.max_abs_error_ms) else np.nan
            ),
            "tolerance_ms": round(self.tolerance_ms, 2),
            "sampling_interval_ms": round(self.sampling_interval_ms, 4),
        }


def match_beats_one_to_one(
    ref_peaks: np.ndarray | list[int],
    det_peaks: np.ndarray | list[int] | list[float],
    fs: float,
    tolerance_ms: float = 150.0,
) -> BeatEvaluationResult:
    """Compare detected QRS peaks against expert reference annotations (strict 1-to-1).

    Preserves the original sorted reference indices so callers can stratify TP/FN
    by beat morphology symbol (e.g., N, A, V, Q). Supports both discrete integer
    sample indices (default) and floating-point sub-sample peak coordinates (e.g.,
    from 3-point parabolic vertex interpolation).
    """
    if not np.isfinite(fs) or fs <= 0.0:
        raise ValueError(f"Sampling frequency 'fs' must be finite and > 0 Hz; got {fs!r}.")
    if not np.isfinite(tolerance_ms) or tolerance_ms < 0.0:
        raise ValueError(f"tolerance_ms must be finite and >= 0; got {tolerance_ms!r}.")

    ref_raw = np.asarray(ref_peaks, dtype=np.int64)
    ref_order = np.argsort(ref_raw, kind="mergesort")
    ref = ref_raw[ref_order]
    # Deduplicate exact sample collisions if any while retaining original index mapping
    if ref.size > 1:
        unique_mask = np.concatenate(([True], np.diff(ref) != 0))
        ref = ref[unique_mask]
        ref_order = ref_order[unique_mask]

    det_raw = np.asarray(det_peaks)
    is_float_det = bool(np.issubdtype(det_raw.dtype, np.floating))
    if is_float_det:
        det = np.unique(det_raw.astype(np.float64))
        tol_limit: float | int = (float(tolerance_ms) / 1000.0) * float(fs)
    else:
        det = np.unique(det_raw.astype(np.int64))
        tol_limit = int(round((tolerance_ms / 1000.0) * fs))
    dt_ms = 1000.0 / float(fs)

    candidates: list[tuple[float, int, int]] = []
    if ref.size > 0 and det.size > 0:
        left_indices = np.searchsorted(det, ref - tol_limit, side="left")
        right_indices = np.searchsorted(det, ref + tol_limit, side="right")
        for r_idx in range(ref.size):
            r_samp = float(ref[r_idx]) if is_float_det else int(ref[r_idx])
            for d_idx in range(int(left_indices[r_idx]), int(right_indices[r_idx])):
                d_val = float(det[d_idx]) if is_float_det else int(det[d_idx])
                diff = abs(d_val - r_samp)
                if diff <= (tol_limit + (1e-9 if is_float_det else 0)):
                    candidates.append((float(diff), r_idx, d_idx))

    candidates.sort(key=lambda item: (item[0], item[1], item[2]))

    ref_used = np.zeros(ref.size, dtype=bool)
    det_used = np.zeros(det.size, dtype=bool)
    matched_pairs: list[tuple[int, int]] = []

    for _, r_idx, d_idx in candidates:
        if not ref_used[r_idx] and not det_used[d_idx]:
            ref_used[r_idx] = True
            det_used[d_idx] = True
            matched_pairs.append((r_idx, d_idx))

    matched_pairs.sort(key=lambda pair: pair[0])

    if matched_pairs:
        m_ref_local = np.array([p[0] for p in matched_pairs], dtype=np.int64)
        m_det_local = np.array([p[1] for p in matched_pairs], dtype=np.int64)
        matched_ref = ref[m_ref_local]
        matched_det = det[m_det_local]
        matched_orig_indices = ref_order[m_ref_local]
        # Signed error convention: e_i = t_det,i - t_ref,i
        timing_errors_ms = ((matched_det - matched_ref).astype(np.float64) / fs) * 1000.0
        abs_errors_ms = np.abs(timing_errors_ms)
        mae_ms = float(np.mean(abs_errors_ms))
        mean_signed_ms = float(np.mean(timing_errors_ms))
        std_ms = float(np.std(timing_errors_ms))
        median_abs_ms = float(np.median(abs_errors_ms))
        max_abs_ms = float(np.max(abs_errors_ms))
    else:
        matched_ref = np.empty(0, dtype=np.int64)
        matched_det = np.empty(0, dtype=np.int64)
        matched_orig_indices = np.empty(0, dtype=np.int64)
        timing_errors_ms = np.empty(0, dtype=np.float64)
        mae_ms = float("nan")
        mean_signed_ms = float("nan")
        std_ms = float("nan")
        median_abs_ms = float("nan")
        max_abs_ms = float("nan")

    unmatched_ref = ref[~ref_used]
    unmatched_orig_indices = ref_order[~ref_used]
    unmatched_det = det[~det_used]

    tp = int(matched_ref.size)
    fp = int(unmatched_det.size)
    fn = int(unmatched_ref.size)

    sensitivity = float(tp / (tp + fn)) if (tp + fn) > 0 else 0.0
    ppv = float(tp / (tp + fp)) if (tp + fp) > 0 else 0.0
    f1_score = (
        float(2.0 * sensitivity * ppv / (sensitivity + ppv))
        if (sensitivity + ppv) > 0.0
        else 0.0
    )

    return BeatEvaluationResult(
        n_ref=int(ref.size),
        n_det=int(det.size),
        tp=tp,
        fp=fp,
        fn=fn,
        sensitivity=sensitivity,
        ppv=ppv,
        f1_score=f1_score,
        mae_ms=mae_ms,
        mean_signed_error_ms=mean_signed_ms,
        std_error_ms=std_ms,
        median_abs_error_ms=median_abs_ms,
        max_abs_error_ms=max_abs_ms,
        tolerance_ms=float(tolerance_ms),
        sampling_interval_ms=dt_ms,
        matched_ref_samples=matched_ref,
        matched_det_samples=matched_det,
        matched_ref_indices=matched_orig_indices,
        timing_errors_ms=timing_errors_ms,
        unmatched_ref_samples=unmatched_ref,
        unmatched_ref_indices=unmatched_orig_indices,
        unmatched_det_samples=unmatched_det,
    )


def evaluate_tolerance_sweep(
    ref_peaks: np.ndarray | list[int],
    det_peaks: np.ndarray | list[int],
    fs: float,
    tolerances_ms: tuple[float, ...] | list[float] = DEFAULT_TOLERANCES_MS,
) -> list[dict[str, Any]]:
    """Evaluate strict 1-to-1 detection performance across multiple temporal tolerances.

    Addresses Experiment 1 (sensitivity of Se, PPV, and F1 to the matching window).
    """
    rows: list[dict[str, Any]] = []
    for tol_ms in tolerances_ms:
        tol_samples = int(round((float(tol_ms) / 1000.0) * fs))
        res = match_beats_one_to_one(ref_peaks=ref_peaks, det_peaks=det_peaks, fs=fs, tolerance_ms=float(tol_ms))
        rows.append(
            {
                "tolerance_ms": round(float(tol_ms), 2),
                "tolerance_samples": tol_samples,
                **res.to_metrics_dict(),
            }
        )
    return rows


def compute_micro_aggregated_metrics(metrics_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate TP, FP, FN, and N_ref counts across records by scenario (micro-averaging).

    Prevents confusing unweighted averages of per-record percentages with true
    pooled count-based Sensitivity, PPV, and F1-score.
    """
    grouped = (
        metrics_df.groupby(["scenario_key", "scenario", "condition"], sort=False, as_index=False)[
            ["n_ref", "n_det", "tp", "fp", "fn"]
        ]
        .sum()
    )
    se = np.where((grouped["tp"] + grouped["fn"]) > 0, grouped["tp"] / (grouped["tp"] + grouped["fn"]), 0.0)
    ppv = np.where((grouped["tp"] + grouped["fp"]) > 0, grouped["tp"] / (grouped["tp"] + grouped["fp"]), 0.0)
    f1 = np.where((se + ppv) > 0, 2.0 * se * ppv / (se + ppv), 0.0)

    grouped["sensitivity"] = np.round(se, 6)
    grouped["ppv"] = np.round(ppv, 6)
    grouped["f1_score"] = np.round(f1, 6)
    return grouped
