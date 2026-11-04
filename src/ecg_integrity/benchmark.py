"""Independent QRS detector benchmark against WFDB's XQRS algorithm.

Implements Experiment 4:
Evaluates our Pan–Tompkins-inspired implementation side-by-side with WFDB's
``wfdb.processing.xqrs_detect`` on the exact same input signals, ground-truth
reference beat annotations, and temporal matching tolerances.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from wfdb import processing

from ecg_integrity.detector import (
    QRSDetectorConfig,
    QRSDetectorV2Config,
    detect_qrs_pan_tompkins,
    detect_qrs_pan_tompkins_v2,
)
from ecg_integrity.evaluation import DEFAULT_TOLERANCES_MS, evaluate_tolerance_sweep, match_beats_one_to_one


def detect_qrs_wfdb_xqrs(signal: np.ndarray, fs: float, learn: bool = True) -> np.ndarray:
    """Run WFDB's reference XQRS detector quietly and return sorted integer peak indices."""
    sig_arr = np.ascontiguousarray(signal, dtype=np.float64)
    peaks = processing.xqrs_detect(sig=sig_arr, fs=float(fs), learn=learn, verbose=False)
    return np.asarray(peaks, dtype=np.int64)


def compare_detectors_on_scenarios(
    record_id: str,
    scenarios: list[dict[str, Any]],
    r_peaks_ref: np.ndarray,
    fs: float,
    primary_tolerance_ms: float = 150.0,
    tolerances_ms: tuple[float, ...] = DEFAULT_TOLERANCES_MS,
    config: QRSDetectorConfig | None = None,
    include_v2: bool = False,
    v2_config: QRSDetectorV2Config | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare Custom Pan–Tompkins (v1 and optional v2) vs. WFDB XQRS across scenarios and tolerance sweeps."""
    summary_rows: list[dict[str, Any]] = []
    sweep_rows: list[dict[str, Any]] = []

    for scen in scenarios:
        sig = scen["signal"]
        pt_v1_peaks = detect_qrs_pan_tompkins(sig, fs=fs, config=config).r_peaks
        xqrs_peaks = detect_qrs_wfdb_xqrs(sig, fs=fs, learn=True)

        detectors_list: list[tuple[str, np.ndarray]] = [("Pan-Tompkins (Custom)", pt_v1_peaks)]
        if include_v2:
            pt_v2_peaks = detect_qrs_pan_tompkins_v2(sig, fs=fs, config=v2_config).r_peaks
            detectors_list = [
                ("Pan-Tompkins v1", pt_v1_peaks),
                ("Pan-Tompkins v2", pt_v2_peaks),
                ("WFDB XQRS", xqrs_peaks),
            ]
        else:
            detectors_list.append(("WFDB XQRS", xqrs_peaks))

        for det_name, peaks in detectors_list:
            primary_res = match_beats_one_to_one(
                ref_peaks=r_peaks_ref,
                det_peaks=peaks,
                fs=fs,
                tolerance_ms=primary_tolerance_ms,
            )
            summary_rows.append(
                {
                    "record": str(record_id),
                    "detector": det_name,
                    "scenario_key": scen["key"],
                    "scenario": scen["scenario"],
                    "condition": scen["condition"],
                    **primary_res.to_metrics_dict(),
                }
            )

            for sw in evaluate_tolerance_sweep(
                ref_peaks=r_peaks_ref,
                det_peaks=peaks,
                fs=fs,
                tolerances_ms=tolerances_ms,
            ):
                sweep_rows.append(
                    {
                        "record": str(record_id),
                        "detector": det_name,
                        "scenario_key": scen["key"],
                        "scenario": scen["scenario"],
                        **sw,
                    }
                )

    return pd.DataFrame(summary_rows), pd.DataFrame(sweep_rows)
