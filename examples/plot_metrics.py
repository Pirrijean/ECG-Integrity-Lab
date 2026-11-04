"""Generate summary comparison figure from a saved metrics.csv file.

Usage:
    python examples/plot_metrics.py
    python examples/plot_metrics.py --csv results/record_100_ch0/metrics.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from ecg_integrity.reporting import save_figure3_performance_metrics


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot Sensitivity, PPV, F1-score, and MAE_t from metrics.csv")
    parser.add_argument(
        "--csv",
        type=str,
        default="results/record_100_ch0/metrics.csv",
        help="Path to metrics.csv (default: results/record_100_ch0/metrics.csv)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output PNG path (default: metrics_summary.png alongside the CSV)",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Metrics file not found at {csv_path}. Run 'python -m ecg_integrity.cli run --record 100 --channel 0 --seconds 300' first."
        )

    df = pd.read_csv(csv_path)
    out_path = Path(args.output) if args.output else csv_path.parent / "metrics_summary.png"
    rec = str(df["record"].iloc[0]) if "record" in df.columns else "Unknown"
    ch_name = str(df["channel_name"].iloc[0]) if "channel_name" in df.columns else "ch0"

    saved = save_figure3_performance_metrics(
        metrics_df=df,
        output_path=out_path,
        title_prefix=f"QRS Detection Performance Summary (Record {rec}, {ch_name})",
    )
    print(f"Saved metrics summary plot to: {saved}")


if __name__ == "__main__":
    main()
