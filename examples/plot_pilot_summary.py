"""Generate cross-record comparative figures for the 5-record MIT-BIH pilot study.

Usage:
    python examples/plot_pilot_summary.py --csv results/pilot/pilot_metrics.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot multi-record pilot F1-score and MAE_t comparison")
    parser.add_argument(
        "--csv",
        type=str,
        default="results/pilot/pilot_metrics.csv",
        help="Path to aggregated pilot_metrics.csv",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise FileNotFoundError(f"Aggregated pilot CSV not found at {csv_path}.")

    df = pd.read_csv(csv_path, dtype={"record": str})
    records = sorted(df["record"].unique().tolist())
    scenarios_order = [
        "Baseline",
        "Baseline wander",
        "Mains interference",
        "White noise (20 dB)",
        "White noise (10 dB)",
        "White noise (5 dB)",
    ]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8.2), constrained_layout=True)
    x = np.arange(len(scenarios_order))
    colors = ["#0f172a", "#2563eb", "#16a34a", "#d97706", "#dc2626", "#7c3aed"]

    for idx, rec in enumerate(records):
        sub = df[df["record"] == rec].set_index("scenario").reindex(scenarios_order)
        c = colors[idx % len(colors)]
        ax1.plot(
            x,
            sub["f1_score"].to_numpy(dtype=float) * 100.0,
            marker="o",
            linewidth=1.8,
            color=c,
            label=f"Record {rec}",
        )
        ax2.plot(
            x,
            sub["mae_ms"].to_numpy(dtype=float),
            marker="s",
            linewidth=1.8,
            color=c,
            label=f"Record {rec}",
        )

    ax1.set_xticks(x)
    ax1.set_xticklabels(scenarios_order, rotation=14, ha="right")
    ax1.set_ylabel("F1-Score (%)")
    ax1.set_ylim(max(0.0, float(df["f1_score"].min() * 100.0) - 5.0), 101.0)
    ax1.set_title("Multi-Record Pilot — F1-Score Across Controlled Perturbations (30-min Full Records)", fontweight="bold")
    ax1.grid(True, linestyle="--", alpha=0.35)
    ax1.legend(loc="lower left", ncol=len(records), fontsize=9)

    ax2.set_xticks(x)
    ax2.set_xticklabels(scenarios_order, rotation=14, ha="right")
    ax2.set_ylabel("Temporal MAE_t (ms)")
    ax2.axhline(1000.0 / 360.0, color="#64748b", linestyle=":", linewidth=1.0, label="1 Sample Δt (2.78 ms)")
    ax2.set_title("Multi-Record Pilot — Mean Absolute Temporal Error (MAE_t)", fontweight="bold")
    ax2.grid(True, linestyle="--", alpha=0.35)
    ax2.legend(loc="upper left", ncol=len(records) + 1, fontsize=8.5)

    out_fig = csv_path.parent / "pilot_comparison.png"
    fig.savefig(out_fig, dpi=200)
    plt.close(fig)
    print(f"Saved multi-record pilot comparison figure to: {out_fig}")


if __name__ == "__main__":
    main()
