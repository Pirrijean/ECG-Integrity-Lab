"""Command-line interface for ECG-Integrity-Lab experimental protocol execution."""

from __future__ import annotations

import argparse
import platform
import sys
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
import pandas as pd
import scipy
import wfdb

from ecg_integrity import __version__
from ecg_integrity.dataset import load_mitbih_record
from ecg_integrity.detector import QRSDetectorConfig, detect_qrs_pan_tompkins
from ecg_integrity.digital_phantom import evaluate_digital_phantom, generate_1mv_60bpm_phantom
from ecg_integrity.evaluation import match_beats_one_to_one
from ecg_integrity.reporting import (
    save_figure1_conditioning,
    save_figure2_perturbations,
    save_figure3_performance_metrics,
    save_figure4_timing_error_distribution,
    save_run_metadata_json,
    save_scenario_diagnostic_plot,
)
from ecg_integrity.signal_processing import build_protocol_scenarios


def resolve_record_output_dir(output_arg: str | None, record_id: str, channel: int) -> tuple[Path, Path | None]:
    """Resolve per-record output directory and optional multi-record parent summary path.

    - Default (no --output): ``results/record_<record>_ch<channel>``
    - If ``--output results/pilot`` is passed across multiple records, stores each record
      under ``results/pilot/record_<record>_ch<channel>`` AND updates ``results/pilot/pilot_metrics.csv``
      unless ``--output`` already ends with ``record_<record>_ch<channel>``.
    """
    rec_folder_name = f"record_{record_id}_ch{channel}"
    if output_arg is None:
        return Path("results") / rec_folder_name, None

    provided = Path(output_arg)
    if provided.name == rec_folder_name:
        return provided, None
    return provided / rec_folder_name, provided / "pilot_metrics.csv"


def run_protocol(
    record_id: str,
    channel: int = 0,
    seconds: float | None = 300.0,
    output_arg: str | None = None,
    seed: int = 42,
    tolerance_ms: float = 150.0,
    cache_dir: str | None = "data/raw/mitdb",
) -> dict[str, Any]:
    """Execute the end-to-end ECG conditioning, perturbation, detection, and QA protocol."""
    out_dir, aggregate_csv_path = resolve_record_output_dir(output_arg, str(record_id), channel)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Software calibration phantom pre-flight check
    phantom_sig = generate_1mv_60bpm_phantom(duration_s=10.0, fs=360.0, amplitude_mv=1.0, bpm=60.0)
    phantom_check = evaluate_digital_phantom(phantom_sig, fs=360.0, expected_amplitude_mv=1.0, expected_bpm=60.0)

    # 2. Load MIT-BIH record and validate Step 9 metadata
    record_data = load_mitbih_record(
        record_id=record_id,
        channel=channel,
        seconds=seconds,
        cache_dir=cache_dir,
        expected_fs=360.0,
    )

    detector_cfg = QRSDetectorConfig()
    scenarios = build_protocol_scenarios(record_data.signal, fs=record_data.fs, seed=seed)

    rows: list[dict[str, Any]] = []
    signed_errors_list: list[np.ndarray] = []
    scenario_labels: list[str] = []
    baseline_det_output = None

    for scen in scenarios:
        det_out = detect_qrs_pan_tompkins(scen["signal"], fs=record_data.fs, config=detector_cfg)
        if scen["key"] == "baseline":
            baseline_det_output = det_out

        eval_res = match_beats_one_to_one(
            ref_peaks=record_data.r_peaks_ref,
            det_peaks=det_out.r_peaks,
            fs=record_data.fs,
            tolerance_ms=tolerance_ms,
        )

        # Save per-scenario diagnostic plot (baseline.png, baseline_wander.png, ...)
        save_scenario_diagnostic_plot(
            output_path=out_dir / f"{scen['key']}.png",
            record_id=record_data.record_id,
            channel_name=record_data.channel_name,
            scenario_title=scen["scenario"],
            condition_desc=scen["condition"],
            input_signal=scen["signal"],
            det_output=det_out,
            eval_result=eval_res,
            r_peaks_ref=record_data.r_peaks_ref,
            fs=record_data.fs,
        )

        row = {
            "record": record_data.record_id,
            "channel": record_data.channel,
            "channel_name": record_data.channel_name,
            "duration_s": record_data.metadata["duration_seconds"],
            "scenario_key": scen["key"],
            "scenario": scen["scenario"],
            "condition": scen["condition"],
            **eval_res.to_metrics_dict(),
        }
        rows.append(row)
        scenario_labels.append(scen["scenario"])
        signed_errors_list.append(eval_res.timing_errors_ms)

    metrics_df = pd.DataFrame(rows)
    metrics_csv_path = out_dir / "metrics.csv"
    metrics_df.to_csv(metrics_csv_path, index=False)

    # Update aggregated multi-record CSV if running under a shared parent output directory (e.g. results/pilot)
    if aggregate_csv_path is not None:
        aggregate_csv_path.parent.mkdir(parents=True, exist_ok=True)
        if aggregate_csv_path.exists():
            existing_df = pd.read_csv(aggregate_csv_path, dtype={"record": str})
            existing_df = existing_df[
                ~(
                    (existing_df["record"].astype(str) == str(record_data.record_id))
                    & (existing_df["channel"].astype(int) == int(record_data.channel))
                )
            ]
            combined_df = pd.concat([existing_df, metrics_df], ignore_index=True)
        else:
            combined_df = metrics_df
        combined_df.to_csv(aggregate_csv_path, index=False)

    # Save Section 12 Deliverable Figures 1–4
    assert baseline_det_output is not None
    save_figure1_conditioning(
        output_path=out_dir / "fig1_conditioning.png",
        record_id=record_data.record_id,
        channel_name=record_data.channel_name,
        clean_signal=record_data.signal,
        det_output=baseline_det_output,
        r_peaks_ref=record_data.r_peaks_ref,
        fs=record_data.fs,
    )
    save_figure2_perturbations(
        output_path=out_dir / "fig2_perturbations.png",
        record_id=record_data.record_id,
        channel_name=record_data.channel_name,
        scenarios=scenarios,
        r_peaks_ref=record_data.r_peaks_ref,
        fs=record_data.fs,
    )
    save_figure3_performance_metrics(
        metrics_df=metrics_df,
        output_path=out_dir / "fig3_performance.png",
        title_prefix=(
            f"Figure 3 — QRS Detection Performance vs Interference "
            f"(Record {record_data.record_id}, {record_data.channel_name})"
        ),
    )
    save_figure3_performance_metrics(
        metrics_df=metrics_df,
        output_path=out_dir / "metrics_summary.png",
        title_prefix=(
            f"QRS Detection Summary Metrics "
            f"(Record {record_data.record_id}, {record_data.channel_name})"
        ),
    )
    save_figure4_timing_error_distribution(
        output_path=out_dir / "fig4_timing_error.png",
        record_id=record_data.record_id,
        channel_name=record_data.channel_name,
        scenario_labels=scenario_labels,
        signed_errors_by_scenario=signed_errors_list,
        fs=record_data.fs,
    )

    # Persist complete QA & metrology traceability metadata
    metadata_payload: dict[str, Any] = {
        "framework": "ECG-Integrity-Lab",
        "version": __version__,
        "environment": {
            "python_version": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy_version": np.__version__,
            "scipy_version": scipy.__version__,
            "pandas_version": pd.__version__,
            "matplotlib_version": matplotlib.__version__,
            "wfdb_version": wfdb.__version__,
        },
        "software_verification_phantom": phantom_check,
        "record_qa_metadata": record_data.metadata,
        "processing_parameters": {
            "filter_type": "Butterworth bandpass (zero-phase second-order sections via sosfiltfilt)",
            "filter_scope_warning": (
                "5-20 Hz QRS-detection filter only; attenuates P/T waves and distorts ST morphology. "
                "Not for clinical or diagnostic ECG interpretation."
            ),
            **detector_cfg.to_dict(),
        },
        "evaluation_protocol": {
            "matching_algorithm": "Greedy closest-pair strict one-to-one beat matching",
            "tolerance_ms": tolerance_ms,
            "tolerance_samples": int(round((tolerance_ms / 1000.0) * record_data.fs)),
            "random_seed": seed,
            "scenarios": [
                {"key": s["key"], "scenario": s["scenario"], "condition": s["condition"], "parameters": s["parameters"]}
                for s in scenarios
            ],
        },
        "output_directory": str(out_dir.resolve()),
    }
    metadata_json_path = save_run_metadata_json(out_dir / "run_metadata.json", metadata_payload)

    return {
        "output_dir": out_dir,
        "metrics_csv": metrics_csv_path,
        "metadata_json": metadata_json_path,
        "metrics_df": metrics_df,
        "record_metadata": record_data.metadata,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ecg_integrity",
        description="ECG-Integrity-Lab: Open ECG Signal Conditioning & QRS Detector Verification Framework",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_p = subparsers.add_parser("run", help="Run the ECG conditioning & perturbation verification protocol")
    run_p.add_argument("--record", type=str, default="100", help="MIT-BIH record ID (default: 100)")
    run_p.add_argument("--channel", type=int, default=0, help="Zero-based ECG channel index (default: 0)")
    run_p.add_argument(
        "--seconds",
        type=float,
        default=300.0,
        help="Signal duration to analyze in seconds (default: 300; use 1800 for ~30-min full record)",
    )
    run_p.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output base directory (default: results/record_<ID>_ch<CH>)",
    )
    run_p.add_argument("--seed", type=int, default=42, help="Random seed for reproducible white noise (default: 42)")
    run_p.add_argument(
        "--tolerance-ms",
        type=float,
        default=150.0,
        help="One-to-one beat matching tolerance in milliseconds (default: 150.0)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "run":
        result = run_protocol(
            record_id=args.record,
            channel=args.channel,
            seconds=args.seconds,
            output_arg=args.output,
            seed=args.seed,
            tolerance_ms=args.tolerance_ms,
        )
        rec_meta = result["record_metadata"]
        print(
            f"[QA OK] Record {rec_meta['record_id']} | Channel {rec_meta['selected_channel']} "
            f"({rec_meta['channel_name']}) | fs={rec_meta['fs_hz']} Hz | "
            f"Units={rec_meta['physical_unit']} | ADC Gain={rec_meta['adc_gain_counts_per_mv']} counts/mV | "
            f"ADC Baseline={rec_meta['adc_baseline_counts']} | Ref Beats={rec_meta['reference_beat_count']}"
        )
        print(result["metrics_df"][["scenario", "tp", "fp", "fn", "sensitivity", "ppv", "f1_score", "mae_ms"]].to_string(index=False))
        print(f"\nSaved outputs to: {result['output_dir']}")
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
