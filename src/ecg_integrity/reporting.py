"""Reproducible figure generation, CSV metrics export, and JSON QA metadata logging."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ecg_integrity.detector import QRSDetectionOutput
from ecg_integrity.evaluation import BeatEvaluationResult


def _style_ecg_axes(ax: plt.Axes) -> None:
    """Apply clean engineering oscilloscope-like grid styling to a Matplotlib Axes."""
    ax.grid(True, which="major", linestyle="-", linewidth=0.55, alpha=0.35, color="#94a3b8")
    ax.grid(True, which="minor", linestyle=":", linewidth=0.35, alpha=0.22, color="#cbd5e1")
    ax.minorticks_on()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_scenario_diagnostic_plot(
    output_path: str | Path,
    record_id: str,
    channel_name: str,
    scenario_title: str,
    condition_desc: str,
    input_signal: np.ndarray,
    det_output: QRSDetectionOutput,
    eval_result: BeatEvaluationResult,
    r_peaks_ref: np.ndarray,
    fs: float,
    window_seconds: float = 10.0,
) -> Path:
    """Save a two-panel diagnostic plot for one interference condition (e.g. baseline.png)."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    n_plot = min(input_signal.size, int(round(window_seconds * fs)))
    t = np.arange(n_plot, dtype=np.float64) / fs

    ref_win = r_peaks_ref[r_peaks_ref < n_plot]
    tp_win = eval_result.matched_det_samples[eval_result.matched_det_samples < n_plot]
    fp_win = eval_result.unmatched_det_samples[eval_result.unmatched_det_samples < n_plot]
    fn_win = eval_result.unmatched_ref_samples[eval_result.unmatched_ref_samples < n_plot]

    fig, axes = plt.subplots(2, 1, figsize=(12, 6.2), sharex=True, constrained_layout=True)

    # Top panel: Input ECG (clean or perturbed) + Reference beat annotations
    ax0 = axes[0]
    ax0.plot(t, input_signal[:n_plot], color="#0f172a", linewidth=0.95, label="Input ECG (mV)")
    if ref_win.size > 0:
        ax0.scatter(
            ref_win / fs,
            input_signal[ref_win],
            color="#2563eb",
            marker="v",
            s=42,
            zorder=4,
            label=f"Expert Reference (.atr) [n={eval_result.n_ref}]",
        )
    ax0.set_ylabel("Amplitude (mV)")
    ax0.set_title(
        f"MIT-BIH Record {record_id} ({channel_name}) — {scenario_title} [{condition_desc}]",
        fontsize=11,
        fontweight="bold",
    )
    ax0.legend(loc="upper right", frameon=True, fontsize=8.5)
    _style_ecg_axes(ax0)

    # Bottom panel: 5–20 Hz Butterworth conditioned ECG + TP / FP / FN markers
    ax1 = axes[1]
    ax1.plot(
        t,
        det_output.filtered[:n_plot],
        color="#0284c7",
        linewidth=0.95,
        label="Butterworth Bandpass 5–20 Hz (zero-phase SOS)",
    )
    if tp_win.size > 0:
        ax1.scatter(
            tp_win / fs,
            det_output.filtered[tp_win],
            color="#16a34a",
            marker="o",
            s=38,
            zorder=5,
            label=f"True Positive (TP={eval_result.tp})",
        )
    if fp_win.size > 0:
        ax1.scatter(
            fp_win / fs,
            det_output.filtered[fp_win],
            color="#dc2626",
            marker="x",
            s=52,
            linewidths=1.8,
            zorder=6,
            label=f"False Positive (FP={eval_result.fp})",
        )
    if fn_win.size > 0:
        ax1.scatter(
            fn_win / fs,
            det_output.filtered[fn_win],
            color="#d97706",
            marker="^",
            s=52,
            zorder=6,
            label=f"Missed Reference / FN (FN={eval_result.fn})",
        )

    mae_str = f"{eval_result.mae_ms:.2f} ms" if np.isfinite(eval_result.mae_ms) else "N/A"
    summary_text = (
        f"Se = {eval_result.sensitivity * 100:.2f}%   |   "
        f"PPV = {eval_result.ppv * 100:.2f}%   |   "
        f"F1 = {eval_result.f1_score * 100:.2f}%   |   "
        f"MAE_t = {mae_str}"
    )
    ax1.set_title(summary_text, fontsize=10, color="#1e293b")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Filtered (mV)")
    ax1.legend(loc="upper right", frameon=True, fontsize=8.5, ncol=2)
    _style_ecg_axes(ax1)

    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    return out_file


def save_figure1_conditioning(
    output_path: str | Path,
    record_id: str,
    channel_name: str,
    clean_signal: np.ndarray,
    det_output: QRSDetectionOutput,
    r_peaks_ref: np.ndarray,
    fs: float,
    window_seconds: float = 8.0,
) -> Path:
    """Figure 1: Signal Conditioning & Pan–Tompkins Pipeline vs Reference Annotations."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    n_plot = min(clean_signal.size, int(round(window_seconds * fs)))
    t = np.arange(n_plot, dtype=np.float64) / fs
    ref_win = r_peaks_ref[r_peaks_ref < n_plot]
    det_win = det_output.r_peaks[det_output.r_peaks < n_plot]

    fig, axes = plt.subplots(3, 1, figsize=(12, 7.6), sharex=True, constrained_layout=True)

    # Panel A: Original Calibrated ECG
    axes[0].plot(t, clean_signal[:n_plot], color="#0f172a", linewidth=1.0, label="Original Calibrated ECG x[n]")
    if ref_win.size > 0:
        axes[0].scatter(
            ref_win / fs,
            clean_signal[ref_win],
            color="#2563eb",
            marker="v",
            s=45,
            zorder=4,
            label="Expert Beat Annotation (.atr)",
        )
    axes[0].set_ylabel("Amplitude (mV)")
    axes[0].set_title(
        f"Figure 1 — Signal Conditioning & QRS Enhancement Pipeline (Record {record_id}, {channel_name})",
        fontsize=11.5,
        fontweight="bold",
    )
    axes[0].legend(loc="upper right", fontsize=8.5)
    _style_ecg_axes(axes[0])

    # Panel B: Zero-phase Butterworth 5–20 Hz Filtered Signal
    axes[1].plot(
        t,
        det_output.filtered[:n_plot],
        color="#0284c7",
        linewidth=1.0,
        label="Conditioned ECG y[n] (Butterworth 5–20 Hz, Order 4, sosfiltfilt)",
    )
    if det_win.size > 0:
        axes[1].scatter(
            det_win / fs,
            det_output.filtered[det_win],
            color="#16a34a",
            marker="o",
            s=40,
            zorder=5,
            label="Relocated QRS Detections",
        )
    axes[1].set_ylabel("Filtered (mV)")
    axes[1].set_title(
        "QRS-Oriented Bandpass Output (Note: non-diagnostic bandwidth; attenuates P/T & ST morphology)",
        fontsize=9.5,
        color="#334155",
    )
    axes[1].legend(loc="upper right", fontsize=8.5)
    _style_ecg_axes(axes[1])

    # Panel C: Moving-Window Integrated Energy M[n]
    axes[2].plot(
        t,
        det_output.integrated[:n_plot],
        color="#7c3aed",
        linewidth=1.05,
        label="Moving-Window Integrated Energy M[n] (120 ms window)",
    )
    if det_win.size > 0:
        axes[2].scatter(
            det_win / fs,
            det_output.integrated[det_win],
            color="#16a34a",
            marker="o",
            s=36,
            zorder=5,
            label="Detected QRS Envelope Peaks",
        )
    axes[2].set_xlabel("Time (s)")
    axes[2].set_ylabel("Energy (mV²)")
    axes[2].legend(loc="upper right", fontsize=8.5)
    _style_ecg_axes(axes[2])

    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    return out_file


def save_figure2_perturbations(
    output_path: str | Path,
    record_id: str,
    channel_name: str,
    scenarios: list[dict[str, Any]],
    r_peaks_ref: np.ndarray,
    fs: float,
    window_seconds: float = 6.0,
) -> Path:
    """Figure 2: Multi-panel visual comparison of all 6 controlled ECG signal conditions."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    n_rows = len(scenarios)
    n_plot = min(scenarios[0]["signal"].size, int(round(window_seconds * fs)))
    t = np.arange(n_plot, dtype=np.float64) / fs
    ref_win = r_peaks_ref[r_peaks_ref < n_plot]

    palette = ["#0f172a", "#0369a1", "#4f46e5", "#0d9488", "#d97706", "#dc2626"]
    fig, axes = plt.subplots(n_rows, 1, figsize=(12, 9.8), sharex=True, constrained_layout=True)

    for idx, (scen, ax) in enumerate(zip(scenarios, axes)):
        sig = scen["signal"][:n_plot]
        color = palette[idx % len(palette)]
        ax.plot(t, sig, color=color, linewidth=0.88)
        for r_samp in ref_win:
            ax.axvline(r_samp / fs, color="#64748b", linestyle="--", linewidth=0.65, alpha=0.55)
        ax.set_ylabel("mV", fontsize=8.5)
        ax.set_title(
            f"{idx + 1}. {scen['scenario']} — {scen['condition']}",
            loc="left",
            fontsize=9.5,
            fontweight="bold",
        )
        _style_ecg_axes(ax)

    axes[0].text(
        0.995,
        0.88,
        "Dashed vertical lines: expert .atr beat locations",
        transform=axes[0].transAxes,
        ha="right",
        va="top",
        fontsize=8.2,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "#f8fafc", "edgecolor": "#cbd5e1"},
    )
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle(
        f"Figure 2 — Controlled Electrical & Physiological Perturbations (Record {record_id}, {channel_name})",
        fontsize=12,
        fontweight="bold",
    )

    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    return out_file


def save_figure3_performance_metrics(
    metrics_df: pd.DataFrame,
    output_path: str | Path,
    title_prefix: str = "Figure 3 — QRS Detection Performance Under Controlled Perturbations",
) -> Path:
    """Figure 3: Bar chart of Sensitivity, PPV, F1-score, and Temporal Error (MAE_t)."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    labels = metrics_df["scenario"].astype(str).tolist()
    x = np.arange(len(labels))
    width = 0.25

    se_pct = metrics_df["sensitivity"].to_numpy(dtype=np.float64) * 100.0
    ppv_pct = metrics_df["ppv"].to_numpy(dtype=np.float64) * 100.0
    f1_pct = metrics_df["f1_score"].to_numpy(dtype=np.float64) * 100.0
    mae_vals = metrics_df["mae_ms"].to_numpy(dtype=np.float64)

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(11.5, 7.6), gridspec_kw={"height_ratios": [1.45, 1.0]}, constrained_layout=True
    )

    rects1 = ax_top.bar(x - width, se_pct, width, label="Sensitivity (Se %)", color="#2563eb")
    rects2 = ax_top.bar(x, ppv_pct, width, label="Positive Predictive Value (PPV %)", color="#16a34a")
    rects3 = ax_top.bar(x + width, f1_pct, width, label="F1-Score (%)", color="#7c3aed")

    min_metric = float(np.nanmin(np.concatenate([se_pct, ppv_pct, f1_pct])))
    y_floor = max(0.0, min(80.0, np.floor((min_metric - 8.0) / 5.0) * 5.0))
    ax_top.set_ylim(y_floor, 103.5)
    ax_top.set_ylabel("Score (%)")
    ax_top.set_xticks(x)
    ax_top.set_xticklabels(labels, rotation=15, ha="right")
    ax_top.set_title(title_prefix, fontsize=11.5, fontweight="bold")
    ax_top.legend(loc="lower left", frameon=True, fontsize=9)
    _style_ecg_axes(ax_top)

    for rect_group in (rects1, rects2, rects3):
        for rect in rect_group:
            h = rect.get_height()
            ax_top.annotate(
                f"{h:.1f}",
                xy=(rect.get_x() + rect.get_width() / 2.0, h),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7.5,
            )

    # Bottom subplot: Temporal Localization Error MAE_t (ms) + FP/FN annotations
    bars_mae = ax_bot.bar(x, np.nan_to_num(mae_vals, nan=0.0), width=0.45, color="#0284c7", label="MAE_t (ms)")
    ax_bot.axhline(1000.0 / 360.0, color="#64748b", linestyle="--", linewidth=0.9, label="1 Sample Δt (2.78 ms)")
    ax_bot.set_ylabel("Temporal MAE (ms)")
    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels(labels, rotation=15, ha="right")
    max_mae = float(np.nanmax(mae_vals)) if np.any(np.isfinite(mae_vals)) else 10.0
    ax_bot.set_ylim(0.0, max(10.0, max_mae * 1.32))
    ax_bot.set_title("Mean Absolute Temporal Localization Error (MAE_t) & Error Counts (FP / FN)", fontsize=10)
    ax_bot.legend(loc="upper left", frameon=True, fontsize=8.5)
    _style_ecg_axes(ax_bot)

    for idx, bar in enumerate(bars_mae):
        h = bar.get_height()
        fp_c = int(metrics_df.iloc[idx]["fp"])
        fn_c = int(metrics_df.iloc[idx]["fn"])
        ax_bot.annotate(
            f"{h:.2f} ms\n(FP:{fp_c}, FN:{fn_c})",
            xy=(bar.get_x() + bar.get_width() / 2.0, h),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8.0,
        )

    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    return out_file


def save_figure4_timing_error_distribution(
    output_path: str | Path,
    record_id: str,
    channel_name: str,
    scenario_labels: list[str],
    signed_errors_by_scenario: list[np.ndarray],
    fs: float,
) -> Path:
    """Figure 4: Distribution of signed and absolute QRS localization timing error (ms)."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    dt_ms = 1000.0 / fs
    abs_errors_by_scenario = [np.abs(arr) for arr in signed_errors_by_scenario]

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(12.5, 5.4), constrained_layout=True)

    # Left panel: Signed timing error boxplot (t_det - t_ref)
    bp1 = ax_left.boxplot(
        [arr if arr.size > 0 else np.array([0.0]) for arr in signed_errors_by_scenario],
        tick_labels=scenario_labels,
        patch_artist=True,
        showfliers=True,
        flierprops={"marker": ".", "markersize": 3.5, "alpha": 0.45},
    )
    for patch in bp1["boxes"]:
        patch.set_facecolor("#bae6fd")
        patch.set_edgecolor("#0369a1")
    ax_left.axhline(0.0, color="#0f172a", linestyle="-", linewidth=0.9, label="Zero Timing Bias (0 ms)")
    ax_left.axhline(dt_ms, color="#64748b", linestyle=":", linewidth=0.85, label=f"±1 Sample (±{dt_ms:.2f} ms)")
    ax_left.axhline(-dt_ms, color="#64748b", linestyle=":", linewidth=0.85)
    ax_left.set_ylabel("Signed Timing Error t_det − t_ref (ms)")
    ax_left.set_title("Signed Localization Error Distribution", fontsize=10.5, fontweight="bold")
    ax_left.tick_params(axis="x", rotation=22)
    ax_left.legend(loc="upper left", fontsize=8.2)
    _style_ecg_axes(ax_left)

    # Right panel: Absolute timing error |t_det - t_ref| boxplot + mean markers
    bp2 = ax_right.boxplot(
        [arr if arr.size > 0 else np.array([0.0]) for arr in abs_errors_by_scenario],
        tick_labels=scenario_labels,
        patch_artist=True,
        showmeans=True,
        meanprops={"marker": "D", "markerfacecolor": "#dc2626", "markeredgecolor": "#991b1b", "markersize": 5},
        flierprops={"marker": ".", "markersize": 3.5, "alpha": 0.45},
    )
    for patch in bp2["boxes"]:
        patch.set_facecolor("#ddd6fe")
        patch.set_edgecolor("#5b21b6")
    ax_right.axhline(dt_ms, color="#64748b", linestyle="--", linewidth=0.85, label=f"1 Sample Δt ({dt_ms:.2f} ms)")
    ax_right.set_ylabel("Absolute Timing Error |t_det − t_ref| (ms)")
    ax_right.set_title("Absolute Localization Error (Red Diamond = MAE_t)", fontsize=10.5, fontweight="bold")
    ax_right.tick_params(axis="x", rotation=22)
    ax_right.legend(loc="upper left", fontsize=8.2)
    _style_ecg_axes(ax_right)

    fig.suptitle(
        f"Figure 4 — QRS Temporal Localization Error Distribution (Record {record_id}, {channel_name})",
        fontsize=12,
        fontweight="bold",
    )

    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    return out_file


def save_run_metadata_json(
    output_path: str | Path,
    metadata_payload: dict[str, Any],
) -> Path:
    """Persist structured JSON metadata for full experimental and QA traceability."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    enriched = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        **metadata_payload,
    }
    out_file.write_text(json.dumps(enriched, indent=2, ensure_ascii=False), encoding="utf-8")
    return out_file


def save_tolerance_sweep_plot(
    sweep_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "Sensitivity of Detection Metrics to Temporal Matching Tolerance",
) -> Path:
    """Experiment 1 plot: Se, PPV, and F1 vs. matching window (8.33 to 150 ms)."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    scenarios = sweep_df["scenario"].drop_duplicates().tolist()
    palette = ["#0f172a", "#0284c7", "#4f46e5", "#0d9488", "#d97706", "#dc2626"]
    markers = ["o", "s", "^", "D", "v", "P"]

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), sharex=True, constrained_layout=True)
    metrics_info = [
        ("sensitivity", "Sensitivity (%)"),
        ("ppv", "Positive Predictive Value (%)"),
        ("f1_score", "F1-Score (%)"),
    ]

    for col, (col_name, ylabel) in enumerate(metrics_info):
        ax = axes[col]
        for idx, scen in enumerate(scenarios):
            sub = sweep_df[sweep_df["scenario"] == scen].sort_values("tolerance_ms")
            ax.plot(
                sub["tolerance_ms"].to_numpy(dtype=float),
                sub[col_name].to_numpy(dtype=float) * 100.0,
                marker=markers[idx % len(markers)],
                color=palette[idx % len(palette)],
                linewidth=1.7,
                markersize=5,
                label=scen,
            )
        ax.set_xlabel("Matching Tolerance ±Δt (ms)")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel, fontsize=10.5, fontweight="bold")
        ax.set_xticks([8.33, 25, 50, 100, 150])
        _style_ecg_axes(ax)

    axes[0].legend(loc="lower right", fontsize=8.0, frameon=True)
    fig.suptitle(title, fontsize=12, fontweight="bold")
    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    return out_file


def save_nstdb_stress_test_plot(
    nstdb_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "MIT-BIH Noise Stress Test (bw, ma, em) — Real Acquisition Noise",
) -> Path:
    """Experiment 2 plot: F1-score, Se, PPV, and MAE_t across physical noise SNRs (24 to -6 dB)."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    noise_colors = {"bw": "#0284c7", "ma": "#d97706", "em": "#dc2626"}
    noise_markers = {"bw": "o", "ma": "s", "em": "^"}

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 7.6), sharex=True, constrained_layout=True)
    panels = [
        (axes[0, 0], "sensitivity", "Sensitivity (%)", 100.0),
        (axes[0, 1], "ppv", "Positive Predictive Value (%)", 100.0),
        (axes[1, 0], "f1_score", "F1-Score (%)", 100.0),
        (axes[1, 1], "mae_ms", "Temporal MAE_t (ms)", 1.0),
    ]

    for ax, col_name, ylabel, scale in panels:
        for ntype in ("bw", "ma", "em"):
            sub = nstdb_df[nstdb_df["noise_type"] == ntype].sort_values("snr_db", ascending=False)
            if sub.empty:
                continue
            lbl = str(sub["noise_label"].iloc[0])
            ax.plot(
                sub["snr_db"].to_numpy(dtype=float),
                sub[col_name].to_numpy(dtype=float) * scale,
                marker=noise_markers[ntype],
                color=noise_colors[ntype],
                linewidth=1.8,
                label=lbl,
            )
        ax.invert_xaxis()  # High SNR (24 dB) on left -> Severe noise (-6 dB) on right
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel, fontsize=10.5, fontweight="bold")
        _style_ecg_axes(ax)

    axes[1, 0].set_xlabel("Target Signal-to-Noise Ratio (dB) [Decreasing →]")
    axes[1, 1].set_xlabel("Target Signal-to-Noise Ratio (dB) [Decreasing →]")
    axes[1, 1].axhline(1000.0 / 360.0, color="#64748b", linestyle=":", linewidth=1.0, label="1 Sample Δt (2.78 ms)")
    axes[0, 0].legend(loc="lower left", fontsize=8.5, frameon=True)
    axes[1, 1].legend(loc="upper left", fontsize=8.0, frameon=True)
    fig.suptitle(title, fontsize=12, fontweight="bold")
    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    return out_file


def save_pvc_morphology_comparison_plot(
    morph_df: pd.DataFrame,
    output_path: str | Path,
    record_id: str,
    condition_label: str,
) -> Path:
    """Experiment 3 plot: Compare A_QRS, W_QRS, max |dV/dt|, and E_QRS for TP vs FN beats."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    pvc_df = morph_df[morph_df["symbol"] == "V"]
    target_df = pvc_df if (not pvc_df.empty and (pvc_df["outcome"] == "FN").any()) else morph_df
    subset_name = "Ventricular Ectopic Beats (V / PVC)" if target_df is pvc_df else "All Reference Beats"

    tp_df = target_df[target_df["outcome"] == "TP"]
    fn_df = target_df[target_df["outcome"] == "FN"]

    features = [
        ("a_qrs_mv", "A_QRS Peak-to-Peak (mV)"),
        ("w_qrs_ms", "Estimated Width W_QRS (ms)"),
        ("max_dvdt_filt_mv_s", "Max Filtered Slope |dV/dt| (mV/s)"),
        ("e_qrs_over_threshold", "Integrated Energy / Threshold (E_QRS / I1)"),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(13.6, 4.6), constrained_layout=True)
    for ax, (col, label) in zip(axes, features):
        tp_vals = tp_df[col].dropna().to_numpy(dtype=float)
        fn_vals = fn_df[col].dropna().to_numpy(dtype=float)
        data = [tp_vals if tp_vals.size > 0 else np.array([0.0]), fn_vals if fn_vals.size > 0 else np.array([0.0])]

        bp = ax.boxplot(
            data,
            tick_labels=[f"TP (n={tp_vals.size})", f"FN (n={fn_vals.size})"],
            patch_artist=True,
            showmeans=True,
        )
        bp["boxes"][0].set_facecolor("#bbf7d0")
        bp["boxes"][0].set_edgecolor("#15803d")
        bp["boxes"][1].set_facecolor("#fecaca")
        bp["boxes"][1].set_edgecolor("#b91c1c")
        if col == "e_qrs_over_threshold":
            ax.axhline(1.0, color="#dc2626", linestyle="--", linewidth=1.0, label="Threshold I1 = 1.0")
            ax.legend(loc="upper right", fontsize=8.0)
        ax.set_title(label, fontsize=9.5, fontweight="bold")
        _style_ecg_axes(ax)

    fig.suptitle(
        f"Record {record_id} ({condition_label}) — Biophysical Feature Comparison: {subset_name}",
        fontsize=11.5,
        fontweight="bold",
    )
    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    return out_file


def save_xqrs_benchmark_plot(
    bench_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "Detector Benchmark: Custom Pan–Tompkins vs. WFDB XQRS",
) -> Path:
    """Experiment 4 / v0.3.0 plot: Side-by-side F1-score and MAE_t comparison across 2 or 3 detectors."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    scenarios = bench_df["scenario"].drop_duplicates().tolist()
    detectors = bench_df["detector"].drop_duplicates().tolist()
    n_det = len(detectors)
    x = np.arange(len(scenarios))
    width = 0.76 / max(1, n_det)

    color_map = {
        "Pan-Tompkins (Custom)": "#2563eb",
        "Pan-Tompkins v1": "#64748b",
        "Pan-Tompkins v2": "#16a34a",
        "WFDB XQRS": "#0284c7",
    }
    fallback_colors = ["#2563eb", "#16a34a", "#0d9488", "#7c3aed"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12.2, 7.6), sharex=True, constrained_layout=True)

    for idx, det_name in enumerate(detectors):
        offset = (idx - (n_det - 1) / 2.0) * width
        sub = bench_df[bench_df["detector"] == det_name].set_index("scenario").reindex(scenarios)
        c = color_map.get(det_name, fallback_colors[idx % len(fallback_colors)])
        ax1.bar(x + offset, sub["f1_score"].to_numpy(dtype=float) * 100.0, width, label=det_name, color=c)
        ax2.bar(x + offset, sub["mae_ms"].to_numpy(dtype=float), width, label=det_name, color=c)

    min_f1 = float(np.nanmin(bench_df["f1_score"].to_numpy(dtype=float) * 100.0))
    ax1.set_ylim(max(0.0, min(95.0, min_f1 - 2.5)), 100.8)
    ax1.set_ylabel("F1-Score (%)")
    ax1.set_title(title, fontsize=11.5, fontweight="bold")
    ax1.legend(loc="lower left", fontsize=9, ncol=n_det)
    _style_ecg_axes(ax1)

    ax2.axhline(1000.0 / 360.0, color="#dc2626", linestyle="--", linewidth=0.9, label="1 Sample Δt (2.78 ms)")
    ax2.set_ylabel("Temporal MAE_t (ms)")
    ax2.set_xticks(x)
    ax2.set_xticklabels(scenarios, rotation=16, ha="right")
    ax2.set_title("Temporal Localization Accuracy Comparison (MAE_t)", fontsize=10.5, fontweight="bold")
    ax2.legend(loc="upper left", fontsize=8.5, ncol=n_det + 1)
    _style_ecg_axes(ax2)

    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    return out_file


def save_v03_three_way_summary_plot(
    pilot_3way_df: pd.DataFrame,
    output_path: str | Path,
) -> Path:
    """v0.3.0 Overview Plot: 3-way comparison (PT v1 vs PT v2 vs WFDB XQRS) across all 5 pilot records."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    detectors = ["Pan-Tompkins v1", "Pan-Tompkins v2", "WFDB XQRS"]
    colors = {"Pan-Tompkins v1": "#64748b", "Pan-Tompkins v2": "#16a34a", "WFDB XQRS": "#0284c7"}
    records = ["100", "101", "103", "105", "106"]

    base_df = pilot_3way_df[pilot_3way_df["scenario_key"] == "baseline"]
    awgn5_df = pilot_3way_df[pilot_3way_df["scenario_key"] == "white_noise_5dB"]

    fig, axes = plt.subplots(2, 2, figsize=(13.6, 8.8), constrained_layout=True)

    # Panel A: Baseline MAE_t across all 5 records (Target 1: Sub-millisecond R-peak alignment)
    ax = axes[0, 0]
    x = np.arange(len(records))
    w = 0.25
    for idx, det in enumerate(detectors):
        sub = base_df[base_df["detector"] == det].set_index("record").reindex(records)
        bars = ax.bar(x + (idx - 1) * w, sub["mae_ms"].to_numpy(dtype=float), w, label=det, color=colors[det])
        for b in bars:
            h = b.get_height()
            ax.annotate(f"{h:.2f}", (b.get_x() + b.get_width() / 2, h), xytext=(0, 2), textcoords="offset points", ha="center", fontsize=7.2)
    ax.axhline(1000.0 / 360.0, color="#dc2626", linestyle=":", linewidth=0.95, label="1 Sample Δt (2.78 ms)")
    ax.set_xticks(x)
    ax.set_xticklabels([f"Rec {r}" for r in records])
    ax.set_ylabel("Baseline MAE_t (ms)")
    ax.set_ylim(0.0, 3.05)
    ax.set_title("A. Target 1 — R-Peak Timing Error MAE_t (Baseline, 30m)", fontsize=10.5, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8.0, ncol=2)
    _style_ecg_axes(ax)

    # Panel B: Record 105 False Positives across Baseline, AWGN 5 dB, and Real EM 6 dB (Target 2: Motion SQI Gate)
    ax = axes[0, 1]
    r105_keys = ["baseline", "white_noise_5dB", "nstdb_em_6dB"]
    r105_lbls = ["Baseline", "AWGN 5 dB", "Real EM (+6 dB)"]
    r105_df = pilot_3way_df[(pilot_3way_df["record"] == "105") & (pilot_3way_df["scenario_key"].isin(r105_keys))]
    x105 = np.arange(len(r105_keys))
    for idx, det in enumerate(detectors):
        sub = r105_df[r105_df["detector"] == det].set_index("scenario_key").reindex(r105_keys)
        bars = ax.bar(x105 + (idx - 1) * w, sub["fp"].to_numpy(dtype=float), w, label=det, color=colors[det])
        for b in bars:
            h = b.get_height()
            ax.annotate(f"{int(round(h))}", (b.get_x() + b.get_width() / 2, h), xytext=(0, 2), textcoords="offset points", ha="center", fontsize=8.0, fontweight="bold")
    ax.set_xticks(x105)
    ax.set_xticklabels(r105_lbls)
    ax.set_ylabel("False Positives (FP count)")
    ax.set_ylim(0, 52)
    ax.set_title("B. Target 2 — Motion Artifact False Positives (Record 105)", fontsize=10.5, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8.0)
    _style_ecg_axes(ax)

    # Panel C: Record 106 Missed Beats (FN) across Baseline, AWGN 5 dB, and Real EM 6 dB (Target 3: PVC Rescue)
    ax = axes[1, 0]
    r106_df = pilot_3way_df[(pilot_3way_df["record"] == "106") & (pilot_3way_df["scenario_key"].isin(r105_keys))]
    for idx, det in enumerate(detectors):
        sub = r106_df[r106_df["detector"] == det].set_index("scenario_key").reindex(r105_keys)
        bars = ax.bar(x105 + (idx - 1) * w, sub["fn"].to_numpy(dtype=float), w, label=det, color=colors[det])
        for b in bars:
            h = b.get_height()
            ax.annotate(f"{int(round(h))}", (b.get_x() + b.get_width() / 2, h), xytext=(0, 2), textcoords="offset points", ha="center", fontsize=8.0, fontweight="bold")
    ax.set_xticks(x105)
    ax.set_xticklabels(r105_lbls)
    ax.set_ylabel("Missed Beats / False Negatives (FN count)")
    ax.set_ylim(0, 108)
    ax.set_title("C. Target 3 — Missed Multiform PVCs / FNs (Record 106)", fontsize=10.5, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8.0)
    _style_ecg_axes(ax)

    # Panel D: Pooled 5-Record Error Totals (FP + FN out of 10,789 beats) at Baseline & AWGN 5 dB
    ax = axes[1, 1]
    cond_dfs = [("Pooled Baseline\n(10,789 beats)", base_df), ("Pooled AWGN 5 dB\n(10,789 beats)", awgn5_df)]
    xp = np.arange(len(cond_dfs))
    for idx, det in enumerate(detectors):
        tot_errs = []
        annotations = []
        for _, cdf in cond_dfs:
            sub = cdf[cdf["detector"] == det]
            fp_sum = int(sub["fp"].sum())
            fn_sum = int(sub["fn"].sum())
            tot_errs.append(fp_sum + fn_sum)
            annotations.append(f"{fp_sum + fn_sum}\n({fp_sum}FP/{fn_sum}FN)")
        bars = ax.bar(xp + (idx - 1) * w, tot_errs, w, label=det, color=colors[det])
        for b, txt in zip(bars, annotations):
            h = b.get_height()
            ax.annotate(txt, (b.get_x() + b.get_width() / 2, h), xytext=(0, 2), textcoords="offset points", ha="center", fontsize=7.5)
    ax.set_xticks(xp)
    ax.set_xticklabels([c[0] for c in cond_dfs])
    ax.set_ylabel("Total Classification Errors (FP + FN)")
    ax.set_ylim(0, max(165, int(ax.get_ylim()[1] * 1.25)))
    ax.set_title("D. Pooled Total Errors (FP + FN) Across All 5 Pilot Records", fontsize=10.5, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8.0)
    _style_ecg_axes(ax)

    fig.suptitle(
        "ECG-Integrity-Lab v0.3.0 — Three-Way Benchmark: Pan–Tompkins v1 vs. Pan–Tompkins v2 vs. WFDB XQRS",
        fontsize=12.5,
        fontweight="bold",
    )
    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    return out_file


def save_v03_nstdb_comparison_plot(
    nstdb_3way_df: pd.DataFrame,
    output_path: str | Path,
    record_id: str = "100",
) -> Path:
    """v0.3.0 Real Noise Stress Test 3-way comparison (PT v1 vs PT v2 vs WFDB XQRS) on ma and em."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    detectors = ["Pan-Tompkins v1", "Pan-Tompkins v2", "WFDB XQRS"]
    colors = {"Pan-Tompkins v1": "#64748b", "Pan-Tompkins v2": "#16a34a", "WFDB XQRS": "#0284c7"}
    markers = {"Pan-Tompkins v1": "o", "Pan-Tompkins v2": "s", "WFDB XQRS": "^"}

    fig, axes = plt.subplots(2, 2, figsize=(12.6, 7.8), sharex=True, constrained_layout=True)

    for col_idx, (ntype, ntitle) in enumerate([("ma", "Muscle Artifact (nstdb/ma)"), ("em", "Electrode Motion (nstdb/em)")]):
        ax_f1 = axes[0, col_idx]
        ax_mae = axes[1, col_idx]
        sub_n = nstdb_3way_df[nstdb_3way_df["noise_type"] == ntype]

        for det in detectors:
            sub_d = sub_n[sub_n["detector"] == det].sort_values("snr_db", ascending=False)
            snrs = sub_d["snr_db"].to_numpy(dtype=float)
            ax_f1.plot(snrs, sub_d["f1_score"].to_numpy(dtype=float) * 100.0, marker=markers[det], color=colors[det], linewidth=1.85, label=det)
            ax_mae.plot(snrs, sub_d["mae_ms"].to_numpy(dtype=float), marker=markers[det], color=colors[det], linewidth=1.85, label=det)

        ax_f1.set_ylabel("F1-Score (%)")
        ax_f1.set_title(f"{ntitle} — F1-Score (%)", fontsize=10.5, fontweight="bold")
        ax_f1.legend(loc="lower left", fontsize=8.2)
        _style_ecg_axes(ax_f1)

        ax_mae.set_xlim(25.5, -7.5)
        ax_mae.axhline(1000.0 / 360.0, color="#dc2626", linestyle=":", linewidth=0.95, label="1 Sample Δt (2.78 ms)")
        ax_mae.set_xlabel("Target Signal-to-Noise Ratio (dB) [Decreasing →]")
        ax_mae.set_ylabel("Temporal MAE_t (ms)")
        ax_mae.set_title(f"{ntitle} — Temporal MAE_t (ms)", fontsize=10.5, fontweight="bold")
        ax_mae.legend(loc="upper left", fontsize=8.0)
        _style_ecg_axes(ax_mae)

    fig.suptitle(
        f"MIT-BIH Noise Stress Test (Record {record_id}, 30m) — Pan–Tompkins v1 vs. v2 vs. WFDB XQRS",
        fontsize=12,
        fontweight="bold",
    )
    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    return out_file


def save_v03_hrv_and_subsample_plot(
    hrv_df: pd.DataFrame,
    output_path: str | Path,
) -> Path:
    """v0.3.0 Metrology Plot: MAE_t vs quantized .atr and downstream RMSSD HRV bias on Record 100."""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    methods = [
        "Pan-Tompkins v1 (Raw Snap)",
        "Pan-Tompkins v2 (Discrete 5–20 Hz)",
        "Pan-Tompkins v2 (Parabolic Sub-Sample)",
        "WFDB XQRS",
    ]
    colors = {
        "Pan-Tompkins v1 (Raw Snap)": "#64748b",
        "Pan-Tompkins v2 (Discrete 5–20 Hz)": "#16a34a",
        "Pan-Tompkins v2 (Parabolic Sub-Sample)": "#7c3aed",
        "WFDB XQRS": "#0284c7",
    }
    markers = {
        "Pan-Tompkins v1 (Raw Snap)": "o",
        "Pan-Tompkins v2 (Discrete 5–20 Hz)": "s",
        "Pan-Tompkins v2 (Parabolic Sub-Sample)": "D",
        "WFDB XQRS": "^",
    }
    conds = ["Baseline", "White noise (20 dB)", "White noise (10 dB)", "White noise (5 dB)"]
    short_lbls = ["Baseline", "AWGN 20 dB", "AWGN 10 dB", "AWGN 5 dB"]
    x = np.arange(len(conds))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.2), constrained_layout=True)

    for m in methods:
        sub = hrv_df[hrv_df["method"] == m].set_index("scenario").reindex(conds)
        ax1.plot(x, sub["mae_ms"].to_numpy(dtype=float), marker=markers[m], color=colors[m], linewidth=1.9, markersize=6, label=m)
        ax2.plot(
            x,
            np.abs(sub["rmssd_error_percent"].to_numpy(dtype=float)),
            marker=markers[m],
            color=colors[m],
            linewidth=1.9,
            markersize=6,
            label=m,
        )

    ax1.set_xticks(x)
    ax1.set_xticklabels(short_lbls)
    ax1.set_ylabel("MAE_t vs. Discrete .atr Reference (ms)")
    ax1.set_title("A. R-Peak MAE_t Relative to Integer .atr Annotations", fontsize=10.5, fontweight="bold")
    ax1.legend(loc="upper left", fontsize=8.0)
    _style_ecg_axes(ax1)

    ax2.set_xticks(x)
    ax2.set_xticklabels(short_lbls)
    ax2.set_ylabel("|RMSSD Error| Relative to Clean Reference (%)")
    ax2.set_title("B. Downstream Sinus NN-HRV (RMSSD) Magnitude Error (%)", fontsize=10.5, fontweight="bold")
    ax2.legend(loc="upper left", fontsize=8.0)
    _style_ecg_axes(ax2)

    fig.suptitle(
        "Record 100 (n = 2,196 Sinus NN Pairs) — Quantized .atr Alignment vs. Continuous HRV (RMSSD) Fidelity",
        fontsize=11.5,
        fontweight="bold",
    )
    fig.savefig(out_file, dpi=200)
    plt.close(fig)
    return out_file


