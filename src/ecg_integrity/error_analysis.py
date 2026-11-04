"""Morphology-stratified error analysis and False Positive / False Negative diagnostics.

Implements Experiment 3 and Section 3 investigative tools:
1. Stratify Sensitivity and MAE_t by PhysioNet beat symbol (N, A, V, Q, ...).
2. Extract per-beat biophysical/signal features around every reference QRS complex:
   - A_QRS: peak-to-peak amplitude (mV)
   - W_QRS: estimated QRS complex width at half-maximum on |filtered| (ms)
   - max |dV/dt|: peak absolute slope on raw and filtered waveforms (mV/s)
   - E_QRS: peak moving-window integrated energy M[n] (mV^2) and ratio to local threshold
   enabling direct quantitative comparison between detected (TP) and missed (FN) PVCs.
3. Inspect ±300 ms windows around False Positives (e.g., the 33 FPs in Record 105),
   cross-referencing PhysioNet auxiliary rhythm/noise annotations ('~', '|', '+')
   and plotting multi-stage waveforms with the dynamic threshold trace overlaid.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wfdb

from ecg_integrity.dataset import PHYSIONET_BEAT_SYMBOLS
from ecg_integrity.detector import QRSDetectionOutput
from ecg_integrity.evaluation import BeatEvaluationResult


def stratify_by_beat_symbol(
    beat_symbols: list[str],
    eval_result: BeatEvaluationResult,
) -> pd.DataFrame:
    """Break down N_ref, TP, FN, Sensitivity, and MAE_t by PhysioNet beat annotation symbol."""
    matched_err_by_orig_idx: dict[int, float] = {
        int(orig_idx): abs(float(err_ms))
        for orig_idx, err_ms in zip(eval_result.matched_ref_indices, eval_result.timing_errors_ms)
    }

    buckets: dict[str, dict[str, Any]] = {}
    for idx, sym in enumerate(beat_symbols):
        entry = buckets.setdefault(sym, {"symbol": sym, "n_ref": 0, "tp": 0, "fn": 0, "abs_errors": []})
        entry["n_ref"] += 1
        if idx in matched_err_by_orig_idx:
            entry["tp"] += 1
            entry["abs_errors"].append(matched_err_by_orig_idx[idx])
        else:
            entry["fn"] += 1

    rows: list[dict[str, Any]] = []
    for sym, data in sorted(buckets.items(), key=lambda kv: -kv[1]["n_ref"]):
        n_r = data["n_ref"]
        tp_c = data["tp"]
        fn_c = data["fn"]
        errs = data["abs_errors"]
        rows.append(
            {
                "symbol": sym,
                "n_ref": n_r,
                "tp": tp_c,
                "fn": fn_c,
                "sensitivity": round(tp_c / n_r, 6) if n_r > 0 else 0.0,
                "mae_ms": round(float(np.mean(errs)), 4) if errs else np.nan,
            }
        )
    return pd.DataFrame(rows)


def compute_reference_beat_morphology_table(
    raw_signal: np.ndarray,
    det_output: QRSDetectionOutput,
    r_peaks_ref: np.ndarray,
    beat_symbols: list[str],
    eval_result: BeatEvaluationResult,
    fs: float,
    half_window_ms: float = 80.0,
) -> pd.DataFrame:
    """Compute A_QRS, W_QRS, max |dV/dt|, and E_QRS for every expert reference beat."""
    half_win = max(2, int(round((half_window_ms / 1000.0) * fs)))
    matched_set = {int(i) for i in eval_result.matched_ref_indices}
    n_total = raw_signal.size

    records: list[dict[str, Any]] = []
    for idx, (samp, sym) in enumerate(zip(r_peaks_ref, beat_symbols)):
        c = int(samp)
        left = max(0, c - half_win)
        right = min(n_total, c + half_win + 1)

        seg_raw = raw_signal[left:right]
        seg_filt = det_output.filtered[left:right]
        seg_int = det_output.integrated[left:right]

        # A_QRS: peak-to-peak amplitude (mV)
        a_qrs_mv = float(np.max(seg_raw) - np.min(seg_raw))
        filt_peak_mv = float(np.max(np.abs(seg_filt)))

        # W_QRS: full-width at half-maximum of |filtered| around local peak (ms)
        abs_filt = np.abs(seg_filt)
        pk_local = int(np.argmax(abs_filt))
        half_height = 0.5 * float(abs_filt[pk_local])
        l_ptr = pk_local
        while l_ptr > 0 and abs_filt[l_ptr - 1] >= half_height:
            l_ptr -= 1
        r_ptr = pk_local
        while r_ptr < (abs_filt.size - 1) and abs_filt[r_ptr + 1] >= half_height:
            r_ptr += 1
        w_qrs_ms = ((r_ptr - l_ptr + 1) / fs) * 1000.0

        # max |dV/dt| in mV/s (on raw and filtered signals)
        max_dvdt_raw = float(np.max(np.abs(np.diff(seg_raw))) * fs) if seg_raw.size > 1 else 0.0
        max_dvdt_filt = float(np.max(np.abs(np.diff(seg_filt))) * fs) if seg_filt.size > 1 else 0.0

        # E_QRS: peak integrated energy in M[n] and local adaptive threshold I1
        e_qrs = float(np.max(seg_int))
        thresh_i1 = float(det_output.threshold_trace[min(c, n_total - 1)])
        energy_to_thresh = (e_qrs / thresh_i1) if thresh_i1 > 0 else np.nan

        records.append(
            {
                "beat_index": idx,
                "sample": c,
                "time_s": round(c / fs, 4),
                "symbol": sym,
                "detected_tp": idx in matched_set,
                "outcome": "TP" if (idx in matched_set) else "FN",
                "a_qrs_mv": round(a_qrs_mv, 4),
                "filt_peak_mv": round(filt_peak_mv, 4),
                "w_qrs_ms": round(w_qrs_ms, 2),
                "max_dvdt_raw_mv_s": round(max_dvdt_raw, 2),
                "max_dvdt_filt_mv_s": round(max_dvdt_filt, 2),
                "e_qrs_integrated": round(e_qrs, 8),
                "threshold_i1": round(thresh_i1, 8),
                "e_qrs_over_threshold": round(float(energy_to_thresh), 4) if np.isfinite(energy_to_thresh) else np.nan,
            }
        )
    return pd.DataFrame(records)


def inspect_false_positive_events(
    record_id: str,
    raw_signal: np.ndarray,
    det_output: QRSDetectionOutput,
    eval_result: BeatEvaluationResult,
    fs: float,
    cache_dir: str | Path = "data/raw/mitdb",
    window_ms: float = 300.0,
) -> pd.DataFrame:
    """Catalog every False Positive with ±300 ms signal metrics and nearby .atr auxiliary markers."""
    fp_samples = eval_result.unmatched_det_samples
    half_win = max(1, int(round((window_ms / 1000.0) * fs)))
    n_total = raw_signal.size

    # Load raw .atr annotations to inspect nearby non-beat artifact/quality markers ('~', '|', '+')
    aux_samples = np.empty(0, dtype=np.int64)
    aux_symbols: list[str] = []
    rec_file = Path(cache_dir) / str(record_id)
    if (Path(cache_dir) / f"{record_id}.atr").exists():
        ann = wfdb.rdann(str(rec_file), "atr")
        non_beats = [
            (int(s), str(sym), str(aux).strip("\x00 ").strip())
            for s, sym, aux in zip(ann.sample, ann.symbol, ann.aux_note)
            if sym not in PHYSIONET_BEAT_SYMBOLS
        ]
        if non_beats:
            aux_samples = np.array([nb[0] for nb in non_beats], dtype=np.int64)
            aux_symbols = [f"{nb[1]}({nb[2]})" if nb[2] else nb[1] for nb in non_beats]

    rows: list[dict[str, Any]] = []
    for fp_idx, c_raw in enumerate(fp_samples):
        c = int(c_raw)
        left = max(0, c - half_win)
        right = min(n_total, c + half_win + 1)

        seg_raw = raw_signal[left:right]
        seg_filt = det_output.filtered[left:right]
        seg_int = det_output.integrated[left:right]
        thresh_val = float(det_output.threshold_trace[min(c, n_total - 1)])

        nearest_aux_sym = "none"
        nearest_aux_dist_s = np.nan
        if aux_samples.size > 0:
            dists = np.abs(aux_samples - c)
            best_k = int(np.argmin(dists))
            nearest_aux_dist_s = round(float(dists[best_k]) / fs, 3)
            if nearest_aux_dist_s <= 5.0:
                nearest_aux_sym = aux_symbols[best_k]

        rows.append(
            {
                "fp_rank": fp_idx + 1,
                "sample": c,
                "time_s": round(c / fs, 3),
                "raw_ptp_mv": round(float(np.max(seg_raw) - np.min(seg_raw)), 4),
                "filtered_peak_mv": round(float(np.max(np.abs(seg_filt))), 4),
                "integrated_peak": round(float(np.max(seg_int)), 8),
                "threshold_i1": round(thresh_val, 8),
                "peak_over_threshold": round(float(np.max(seg_int)) / thresh_val, 3) if thresh_val > 0 else np.nan,
                "nearby_non_beat_marker_5s": nearest_aux_sym,
                "nearest_marker_dist_s": nearest_aux_dist_s,
            }
        )
    return pd.DataFrame(rows)


def save_fp_multistage_inspection_plot(
    output_path: str | Path,
    record_id: str,
    raw_signal: np.ndarray,
    det_output: QRSDetectionOutput,
    r_peaks_ref: np.ndarray,
    fp_samples: np.ndarray,
    fs: float,
    window_ms: float = 300.0,
    max_examples: int = 4,
) -> Path:
    """Plot ±300 ms 4-stage inspection panels (raw, filtered, derivative, integrated + threshold) for FPs."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    n_show = min(int(fp_samples.size), max_examples)
    if n_show == 0:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No False Positives detected", ha="center", va="center")
        fig.savefig(out_file, dpi=200)
        plt.close(fig)
        return out_file

    # Select evenly spaced examples across the recording so different artifact episodes are shown
    chosen_indices = np.linspace(0, fp_samples.size - 1, num=n_show, dtype=int)
    chosen_fps = fp_samples[chosen_indices]
    half_win = max(1, int(round((window_ms / 1000.0) * fs)))

    fig, axes = plt.subplots(4, n_show, figsize=(3.6 * n_show, 8.2), sharex="col", constrained_layout=True)
    if n_show == 1:
        axes = axes[:, np.newaxis]

    for col, fp_samp in enumerate(chosen_fps):
        c = int(fp_samp)
        left = max(0, c - half_win)
        right = min(raw_signal.size, c + half_win + 1)
        rel_ms = ((np.arange(left, right) - c) / fs) * 1000.0

        # Nearby reference beats if any inside ±300 ms
        refs_in_win = r_peaks_ref[(r_peaks_ref >= left) & (r_peaks_ref < right)]

        # Row 0: Original ECG
        axes[0, col].plot(rel_ms, raw_signal[left:right], color="#0f172a", linewidth=1.0)
        axes[0, col].axvline(0.0, color="#dc2626", linestyle="--", linewidth=1.0, label="FP @ 0 ms")
        for r_s in refs_in_win:
            axes[0, col].axvline(((r_s - c) / fs) * 1000.0, color="#2563eb", linestyle=":", linewidth=1.1, label="Ref .atr")
        axes[0, col].set_title(f"FP @ t={c / fs:.2f} s\n(sample {c})", fontsize=9.5, fontweight="bold")
        if col == 0:
            axes[0, col].set_ylabel("Raw ECG\n(mV)")

        # Row 1: Filtered 5–20 Hz
        axes[1, col].plot(rel_ms, det_output.filtered[left:right], color="#0284c7", linewidth=1.0)
        axes[1, col].axvline(0.0, color="#dc2626", linestyle="--", linewidth=1.0)
        if col == 0:
            axes[1, col].set_ylabel("Filtered 5–20 Hz\n(mV)")

        # Row 2: Derivative d[n]
        axes[2, col].plot(rel_ms, det_output.differentiated[left:right], color="#0d9488", linewidth=1.0)
        axes[2, col].axvline(0.0, color="#dc2626", linestyle="--", linewidth=1.0)
        if col == 0:
            axes[2, col].set_ylabel("Derivative d[n]\n(mV/sample)")

        # Row 3: Integrated Energy M[n] + Adaptive Threshold I1
        axes[3, col].plot(rel_ms, det_output.integrated[left:right], color="#7c3aed", linewidth=1.1, label="M[n]")
        axes[3, col].plot(
            rel_ms,
            det_output.threshold_trace[left:right],
            color="#dc2626",
            linestyle="-.",
            linewidth=1.15,
            label="Threshold I1",
        )
        axes[3, col].axvline(0.0, color="#dc2626", linestyle="--", linewidth=1.0)
        axes[3, col].set_xlabel("Offset from FP (ms)")
        if col == 0:
            axes[3, col].set_ylabel("Integrated M[n]\n(mV²)")
            axes[3, col].legend(loc="upper right", fontsize=7.8)

        for r in range(4):
            axes[r, col].grid(True, linestyle=":", alpha=0.35)
            axes[r, col].spines["top"].set_visible(False)
            axes[r, col].spines["right"].set_visible(False)

    fig.suptitle(
        f"Record {record_id} — Multi-Stage Inspection (±{window_ms:.0f} ms) of Representative False Positives",
        fontsize=11.5,
        fontweight="bold",
    )
    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    return out_file
