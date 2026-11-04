"""Quantitative QRS detection evaluator with strict one-to-one beat matching.

Implements Steps 13 and 14:
- One-to-one temporal matching within a declared tolerance window (default 150 ms).
- Sensitivity (Se), Positive Predictive Value (PPV), F1-score, and Mean Absolute
  Temporal Error (MAE_t in ms).

Note on standards context:
ANSI/AAMI EC57:2012 (R2020) describes methods for evaluating and reporting cardiac
rhythm algorithms. This module implements a documented one-to-one beat matching
protocol for research/QA verification and does NOT claim full conformity to EC57.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


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
    matched_ref_samples: np.ndarray
    matched_det_samples: np.ndarray
    timing_errors_ms: np.ndarray
    unmatched_ref_samples: np.ndarray
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
        }


def match_beats_one_to_one(
    ref_peaks: np.ndarray | list[int],
    det_peaks: np.ndarray | list[int],
    fs: float,
    tolerance_ms: float = 150.0,
) -> BeatEvaluationResult:
    """Compare detected QRS peaks against expert reference annotations (1-to-1).

    Candidate (reference, detection) pairs within ``tolerance_ms`` are ranked by
    increasing absolute temporal distance ``|t_det - t_ref|`` and assigned greedily
    so that each reference beat and each detected peak can participate in at most
    one True Positive match.
    """
    if not np.isfinite(fs) or fs <= 0.0:
        raise ValueError(f"Sampling frequency 'fs' must be finite and > 0 Hz; got {fs!r}.")
    if not np.isfinite(tolerance_ms) or tolerance_ms < 0.0:
        raise ValueError(f"tolerance_ms must be finite and >= 0; got {tolerance_ms!r}.")

    ref = np.unique(np.asarray(ref_peaks, dtype=np.int64))
    det = np.unique(np.asarray(det_peaks, dtype=np.int64))

    tol_samples = int(round((tolerance_ms / 1000.0) * fs))

    candidates: list[tuple[int, int, int]] = []
    if ref.size > 0 and det.size > 0:
        left_indices = np.searchsorted(det, ref - tol_samples, side="left")
        right_indices = np.searchsorted(det, ref + tol_samples, side="right")
        for r_idx in range(ref.size):
            r_samp = int(ref[r_idx])
            for d_idx in range(int(left_indices[r_idx]), int(right_indices[r_idx])):
                diff = abs(int(det[d_idx]) - r_samp)
                if diff <= tol_samples:
                    candidates.append((diff, r_idx, d_idx))

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
        m_ref_idx = np.array([p[0] for p in matched_pairs], dtype=np.int64)
        m_det_idx = np.array([p[1] for p in matched_pairs], dtype=np.int64)
        matched_ref = ref[m_ref_idx]
        matched_det = det[m_det_idx]
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
        timing_errors_ms = np.empty(0, dtype=np.float64)
        mae_ms = float("nan")
        mean_signed_ms = float("nan")
        std_ms = float("nan")
        median_abs_ms = float("nan")
        max_abs_ms = float("nan")

    unmatched_ref = ref[~ref_used]
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
        matched_ref_samples=matched_ref,
        matched_det_samples=matched_det,
        timing_errors_ms=timing_errors_ms,
        unmatched_ref_samples=unmatched_ref,
        unmatched_det_samples=unmatched_det,
    )
