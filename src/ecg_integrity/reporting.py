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
