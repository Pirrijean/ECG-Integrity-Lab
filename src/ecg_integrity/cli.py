"""Command-line interface for ECG-Integrity-Lab (Phase 1 & Phase 2 protocols)."""

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
from ecg_integrity.benchmark import compare_detectors_on_scenarios
from ecg_integrity.dataset import load_mitbih_record
from ecg_integrity.detector import QRSDetectorConfig, detect_qrs_pan_tompkins
from ecg_integrity.digital_phantom import evaluate_digital_phantom, generate_1mv_60bpm_phantom
from ecg_integrity.error_analysis import (
    compute_reference_beat_morphology_table,
    inspect_false_positive_events,
    save_fp_multistage_inspection_plot,
    stratify_by_beat_symbol,
)
from ecg_integrity.evaluation import DEFAULT_TOLERANCES_MS, evaluate_tolerance_sweep, match_beats_one_to_one
from ecg_integrity.hrv import evaluate_nn_hrv_jitter_impact
from ecg_integrity.nstdb import NSTDB_DEFAULT_SNRS_DB, NSTDB_NOISE_TYPES, load_nstdb_noise, mix_ecg_with_real_noise
from ecg_integrity.reporting import (
    save_figure1_conditioning,
    save_figure2_perturbations,
    save_figure3_performance_metrics,
    save_figure4_timing_error_distribution,
    save_nstdb_stress_test_plot,
    save_pvc_morphology_comparison_plot,
    save_run_metadata_json,
    save_scenario_diagnostic_plot,
    save_tolerance_sweep_plot,
    save_xqrs_benchmark_plot,
)
from ecg_integrity.signal_processing import add_white_noise, build_protocol_scenarios


def resolve_record_output_dir(output_arg: str | None, record_id: str, channel: int) -> tuple[Path, Path | None]:
    """Resolve per-record output directory and optional multi-record parent summary path."""
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
    """Execute the end-to-end ECG conditioning, perturbation, tolerance sweep, and HRV protocol."""
    out_dir, aggregate_csv_path = resolve_record_output_dir(output_arg, str(record_id), channel)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Functional and mathematical software verification via digital reference pulse train
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
    sweep_rows: list[dict[str, Any]] = []
    symbol_rows: list[dict[str, Any]] = []
    hrv_rows: list[dict[str, Any]] = []
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

        # Experiment 1: Temporal tolerance sweep (8.33 to 150 ms)
        for sw in evaluate_tolerance_sweep(
            ref_peaks=record_data.r_peaks_ref,
            det_peaks=det_out.r_peaks,
            fs=record_data.fs,
            tolerances_ms=DEFAULT_TOLERANCES_MS,
        ):
            sweep_rows.append(
                {
                    "record": record_data.record_id,
                    "channel": record_data.channel,
                    "scenario_key": scen["key"],
                    "scenario": scen["scenario"],
                    **sw,
                }
            )

        # Experiment 3: Stratification by PhysioNet beat symbol (N, A, V, Q, ...)
        sym_df = stratify_by_beat_symbol(record_data.beat_symbols, eval_res)
        for _, s_row in sym_df.iterrows():
            symbol_rows.append(
                {
                    "record": record_data.record_id,
                    "scenario_key": scen["key"],
                    "scenario": scen["scenario"],
                    **s_row.to_dict(),
                }
            )

        # Section 7: Time-domain NN HRV (RMSSD) localization-jitter impact
        hrv_res = evaluate_nn_hrv_jitter_impact(
            r_peaks_ref=record_data.r_peaks_ref,
            beat_symbols=record_data.beat_symbols,
            eval_result=eval_res,
            fs=record_data.fs,
        )
        hrv_rows.append(
            {
                "record": record_data.record_id,
                "scenario_key": scen["key"],
                "scenario": scen["scenario"],
                "condition": scen["condition"],
                **hrv_res,
            }
        )

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

    sweep_df = pd.DataFrame(sweep_rows)
    sweep_df.to_csv(out_dir / "tolerance_sweep.csv", index=False)
    save_tolerance_sweep_plot(
        sweep_df=sweep_df,
        output_path=out_dir / "fig_tolerance_sweep.png",
        title=f"Sensitivity to Temporal Matching Window (Record {record_data.record_id}, {record_data.channel_name})",
    )

    pd.DataFrame(symbol_rows).to_csv(out_dir / "morphology_stratified.csv", index=False)
    pd.DataFrame(hrv_rows).to_csv(out_dir / "hrv_jitter_metrics.csv", index=False)

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
        "functional_digital_reference_verification": phantom_check,
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
            "signed_error_convention": "e_i = t_det,i - t_ref,i (positive = detection lags behind reference)",
            "discrete_sampling_interval_ms": round(1000.0 / record_data.fs, 4),
            "primary_tolerance_ms": tolerance_ms,
            "tolerance_sweep_ms": list(DEFAULT_TOLERANCES_MS),
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


def run_nstdb_protocol(
    record_id: str,
    channel: int = 0,
    seconds: float | None = 300.0,
    output_arg: str | None = None,
    tolerance_ms: float = 150.0,
    snrs_db: tuple[float, ...] = NSTDB_DEFAULT_SNRS_DB,
) -> pd.DataFrame:
    """Experiment 2: Evaluate detector against physical noise (bw, ma, em) from PhysioNet nstdb."""
    base_dir = Path(output_arg) if output_arg else Path("results") / "nstdb" / f"record_{record_id}_ch{channel}"
    base_dir.mkdir(parents=True, exist_ok=True)

    rec = load_mitbih_record(record_id=record_id, channel=channel, seconds=seconds)
    cfg = QRSDetectorConfig()

    rows: list[dict[str, Any]] = []
    sweep_rows: list[dict[str, Any]] = []
    hrv_rows: list[dict[str, Any]] = []

    for ntype in NSTDB_NOISE_TYPES:
        noise_seg, noise_meta = load_nstdb_noise(ntype, n_samples=rec.signal.size, channel=0, expected_fs=rec.fs)
        for snr_db in snrs_db:
            noisy_sig, mix_info = mix_ecg_with_real_noise(rec.signal, noise_seg, snr_db=snr_db, fs=rec.fs)
            det_out = detect_qrs_pan_tompkins(noisy_sig, fs=rec.fs, config=cfg)
            eval_res = match_beats_one_to_one(rec.r_peaks_ref, det_out.r_peaks, fs=rec.fs, tolerance_ms=tolerance_ms)

            rows.append(
                {
                    "record": rec.record_id,
                    "noise_type": ntype,
                    "noise_label": noise_meta["label"],
                    "snr_db": float(snr_db),
                    "alpha_scale": mix_info["alpha_scale"],
                    **eval_res.to_metrics_dict(),
                }
            )
            for sw in evaluate_tolerance_sweep(rec.r_peaks_ref, det_out.r_peaks, fs=rec.fs):
                sweep_rows.append(
                    {
                        "record": rec.record_id,
                        "noise_type": ntype,
                        "snr_db": float(snr_db),
                        **sw,
                    }
                )
            hrv_rows.append(
                {
                    "record": rec.record_id,
                    "noise_type": ntype,
                    "snr_db": float(snr_db),
                    **evaluate_nn_hrv_jitter_impact(rec.r_peaks_ref, rec.beat_symbols, eval_res, fs=rec.fs),
                }
            )

    df = pd.DataFrame(rows)
    df.to_csv(base_dir / "nstdb_metrics.csv", index=False)
    pd.DataFrame(sweep_rows).to_csv(base_dir / "nstdb_tolerance_sweep.csv", index=False)
    pd.DataFrame(hrv_rows).to_csv(base_dir / "nstdb_hrv_metrics.csv", index=False)
    save_nstdb_stress_test_plot(
        df,
        output_path=base_dir / "fig_nstdb_stress_test.png",
        title=f"MIT-BIH Noise Stress Test (bw, ma, em) — Record {rec.record_id} ({rec.channel_name})",
    )
    return df


def run_diagnose_protocol(
    record_id: str,
    channel: int = 0,
    seconds: float | None = 1800.0,
    snr_db: float | None = None,
    seed: int = 42,
    tolerance_ms: float = 150.0,
    output_arg: str | None = None,
) -> dict[str, Any]:
    """Experiment 3: Deep FP (±300 ms) inspection and TP-vs-FN biophysical morphology analysis."""
    cond_tag = f"awgn_{snr_db:g}dB" if snr_db is not None else "baseline"
    out_dir = Path(output_arg) if output_arg else Path("results") / "diagnostics" / f"record_{record_id}_{cond_tag}"
    out_dir.mkdir(parents=True, exist_ok=True)

    rec = load_mitbih_record(record_id=record_id, channel=channel, seconds=seconds)
    sig = add_white_noise(rec.signal, snr_db=snr_db, seed=seed, fs=rec.fs) if snr_db is not None else rec.signal.copy()
    cond_label = f"AWGN {snr_db:g} dB" if snr_db is not None else "Baseline"

    det_out = detect_qrs_pan_tompkins(sig, fs=rec.fs)
    eval_res = match_beats_one_to_one(rec.r_peaks_ref, det_out.r_peaks, fs=rec.fs, tolerance_ms=tolerance_ms)

    sym_df = stratify_by_beat_symbol(rec.beat_symbols, eval_res)
    sym_df.to_csv(out_dir / "stratified_by_symbol.csv", index=False)

    morph_df = compute_reference_beat_morphology_table(
        raw_signal=sig,
        det_output=det_out,
        r_peaks_ref=rec.r_peaks_ref,
        beat_symbols=rec.beat_symbols,
        eval_result=eval_res,
        fs=rec.fs,
    )
    morph_df.to_csv(out_dir / "beat_morphology_features.csv", index=False)
    save_pvc_morphology_comparison_plot(morph_df, out_dir / "fig_morphology_tp_vs_fn.png", rec.record_id, cond_label)

    fp_df = inspect_false_positive_events(
        record_id=rec.record_id,
        raw_signal=sig,
        det_output=det_out,
        eval_result=eval_res,
        fs=rec.fs,
    )
    fp_df.to_csv(out_dir / "fp_inspection.csv", index=False)
    save_fp_multistage_inspection_plot(
        output_path=out_dir / "fig_fp_inspection_gallery.png",
        record_id=rec.record_id,
        raw_signal=sig,
        det_output=det_out,
        r_peaks_ref=rec.r_peaks_ref,
        fp_samples=eval_res.unmatched_det_samples,
        fs=rec.fs,
    )

    return {"out_dir": out_dir, "sym_df": sym_df, "morph_df": morph_df, "fp_df": fp_df, "eval": eval_res}


def run_benchmark_protocol(
    record_id: str,
    channel: int = 0,
    seconds: float | None = 300.0,
    seed: int = 42,
    tolerance_ms: float = 150.0,
    output_arg: str | None = None,
) -> pd.DataFrame:
    """Experiment 4: Side-by-side comparison of Custom Pan–Tompkins vs. WFDB XQRS."""
    out_dir = Path(output_arg) if output_arg else Path("results") / "benchmark" / f"record_{record_id}_ch{channel}"
    out_dir.mkdir(parents=True, exist_ok=True)

    rec = load_mitbih_record(record_id=record_id, channel=channel, seconds=seconds)
    scenarios = build_protocol_scenarios(rec.signal, fs=rec.fs, seed=seed)

    # Also append severe physical electrode motion (nstdb/em @ 6 dB) to test real-noise differentiation
    em_noise, _ = load_nstdb_noise("em", n_samples=rec.signal.size, channel=0, expected_fs=rec.fs)
    em_6db, _ = mix_ecg_with_real_noise(rec.signal, em_noise, snr_db=6.0, fs=rec.fs)
    scenarios.append(
        {
            "key": "nstdb_em_6dB",
            "scenario": "Real EM (6 dB)",
            "condition": "NSTDB Electrode Motion 6 dB",
            "signal": em_6db,
            "parameters": {"type": "nstdb_em", "snr_db": 6.0},
        }
    )

    summary_df, sweep_df = compare_detectors_on_scenarios(
        record_id=rec.record_id,
        scenarios=scenarios,
        r_peaks_ref=rec.r_peaks_ref,
        fs=rec.fs,
        primary_tolerance_ms=tolerance_ms,
    )
    summary_df.to_csv(out_dir / "benchmark_xqrs.csv", index=False)
    sweep_df.to_csv(out_dir / "benchmark_tolerance_sweep.csv", index=False)
    save_xqrs_benchmark_plot(
        summary_df,
        output_path=out_dir / "fig_benchmark_xqrs.png",
        title=f"Detector Benchmark: Custom Pan–Tompkins vs. WFDB XQRS (Record {rec.record_id}, {rec.channel_name})",
    )
    return summary_df


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ecg_integrity",
        description="ECG-Integrity-Lab: Open ECG Signal Conditioning & QRS Detector Verification Framework",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand 1: run
    run_p = subparsers.add_parser("run", help="Run synthetic perturbation, tolerance sweep, and HRV protocol")
    run_p.add_argument("--record", type=str, default="100", help="MIT-BIH record ID (default: 100)")
    run_p.add_argument("--channel", type=int, default=0, help="Zero-based ECG channel index (default: 0)")
    run_p.add_argument("--seconds", type=float, default=300.0, help="Duration in seconds (default: 300; 1800 for 30m)")
    run_p.add_argument("--output", type=str, default=None, help="Output directory")
    run_p.add_argument("--seed", type=int, default=42, help="Random seed for AWGN (default: 42)")
    run_p.add_argument("--tolerance-ms", type=float, default=150.0, help="Primary matching tolerance (ms)")

    # Subcommand 2: nstdb (Experiment 2)
    nst_p = subparsers.add_parser("nstdb", help="Run MIT-BIH Noise Stress Test (bw, ma, em across 24..-6 dB SNR)")
    nst_p.add_argument("--record", type=str, default="100", help="MIT-BIH record ID (default: 100)")
    nst_p.add_argument("--channel", type=int, default=0, help="Zero-based ECG channel index")
    nst_p.add_argument("--seconds", type=float, default=300.0, help="Duration in seconds")
    nst_p.add_argument("--output", type=str, default=None, help="Output directory")
    nst_p.add_argument("--tolerance-ms", type=float, default=150.0, help="Primary matching tolerance (ms)")

    # Subcommand 3: diagnose (Experiment 3)
    diag_p = subparsers.add_parser("diagnose", help="Inspect ±300 ms FPs and TP-vs-FN QRS morphology features")
    diag_p.add_argument("--record", type=str, default="105", help="MIT-BIH record ID (e.g. 105 or 106)")
    diag_p.add_argument("--channel", type=int, default=0, help="Zero-based ECG channel index")
    diag_p.add_argument("--seconds", type=float, default=1800.0, help="Duration in seconds (default: 1800)")
    diag_p.add_argument("--snr-db", type=float, default=None, help="Optional AWGN SNR (dB), e.g. 5.0 for Record 106")
    diag_p.add_argument("--output", type=str, default=None, help="Output directory")

    # Subcommand 4: benchmark (Experiment 4)
    bench_p = subparsers.add_parser("benchmark", help="Compare Custom Pan–Tompkins vs. WFDB XQRS detector")
    bench_p.add_argument("--record", type=str, default="100", help="MIT-BIH record ID")
    bench_p.add_argument("--channel", type=int, default=0, help="Zero-based ECG channel index")
    bench_p.add_argument("--seconds", type=float, default=300.0, help="Duration in seconds")
    bench_p.add_argument("--output", type=str, default=None, help="Output directory")

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
            f"Units={rec_meta['physical_unit']} | Ref Beats={rec_meta['reference_beat_count']}"
        )
        print(result["metrics_df"][["scenario", "tp", "fp", "fn", "sensitivity", "ppv", "f1_score", "mae_ms"]].to_string(index=False))
        print(f"Saved outputs to: {result['output_dir']}")
        return 0

    if args.command == "nstdb":
        df = run_nstdb_protocol(
            record_id=args.record,
            channel=args.channel,
            seconds=args.seconds,
            output_arg=args.output,
            tolerance_ms=args.tolerance_ms,
        )
        print(df[["noise_type", "snr_db", "tp", "fp", "fn", "sensitivity", "ppv", "f1_score", "mae_ms"]].to_string(index=False))
        return 0

    if args.command == "diagnose":
        res = run_diagnose_protocol(
            record_id=args.record,
            channel=args.channel,
            seconds=args.seconds,
            snr_db=args.snr_db,
            output_arg=args.output,
        )
        print("=== Beat Symbol Stratification ===")
        print(res["sym_df"].to_string(index=False))
        print(f"False Positives cataloged: {len(res['fp_df'])} | Saved to: {res['out_dir']}")
        return 0

    if args.command == "benchmark":
        df = run_benchmark_protocol(
            record_id=args.record,
            channel=args.channel,
            seconds=args.seconds,
            output_arg=args.output,
        )
        print(df[["detector", "scenario", "tp", "fp", "fn", "sensitivity", "ppv", "f1_score", "mae_ms"]].to_string(index=False))
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
