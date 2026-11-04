"""Public package interface for ECG-Integrity-Lab."""

from __future__ import annotations

__version__ = "0.3.0"

from ecg_integrity.benchmark import compare_detectors_on_scenarios, detect_qrs_wfdb_xqrs
from ecg_integrity.dataset import (
    PHYSIONET_BEAT_SYMBOLS,
    ECGRecordData,
    load_mitbih_record,
    validate_ecg_metadata_and_signal,
)
from ecg_integrity.detector import (
    QRSDetectionOutput,
    QRSDetectionV2Output,
    QRSDetectorConfig,
    QRSDetectorV2Config,
    detect_qrs,
    detect_qrs_pan_tompkins,
    detect_qrs_pan_tompkins_v2,
    refine_peaks_parabolic,
)
from ecg_integrity.digital_phantom import evaluate_digital_phantom, generate_1mv_60bpm_phantom
from ecg_integrity.error_analysis import (
    compute_reference_beat_morphology_table,
    inspect_false_positive_events,
    stratify_by_beat_symbol,
)
from ecg_integrity.evaluation import (
    DEFAULT_TOLERANCES_MS,
    BeatEvaluationResult,
    compute_micro_aggregated_metrics,
    evaluate_tolerance_sweep,
    match_beats_one_to_one,
)
from ecg_integrity.hrv import evaluate_nn_hrv_jitter_impact
from ecg_integrity.nstdb import NSTDB_DEFAULT_SNRS_DB, NSTDB_NOISE_TYPES, load_nstdb_noise, mix_ecg_with_real_noise
from ecg_integrity.signal_processing import (
    add_baseline_wander,
    add_mains_interference,
    add_white_noise,
    bandpass_filter,
    build_protocol_scenarios,
)

__all__ = [
    "__version__",
    "PHYSIONET_BEAT_SYMBOLS",
    "ECGRecordData",
    "load_mitbih_record",
    "validate_ecg_metadata_and_signal",
    "QRSDetectorConfig",
    "QRSDetectorV2Config",
    "QRSDetectionOutput",
    "QRSDetectionV2Output",
    "detect_qrs",
    "detect_qrs_pan_tompkins",
    "detect_qrs_pan_tompkins_v2",
    "refine_peaks_parabolic",
    "generate_1mv_60bpm_phantom",
    "evaluate_digital_phantom",
    "DEFAULT_TOLERANCES_MS",
    "BeatEvaluationResult",
    "match_beats_one_to_one",
    "evaluate_tolerance_sweep",
    "compute_micro_aggregated_metrics",
    "evaluate_nn_hrv_jitter_impact",
    "NSTDB_NOISE_TYPES",
    "NSTDB_DEFAULT_SNRS_DB",
    "load_nstdb_noise",
    "mix_ecg_with_real_noise",
    "stratify_by_beat_symbol",
    "compute_reference_beat_morphology_table",
    "inspect_false_positive_events",
    "detect_qrs_wfdb_xqrs",
    "compare_detectors_on_scenarios",
    "bandpass_filter",
    "add_baseline_wander",
    "add_mains_interference",
    "add_white_noise",
    "build_protocol_scenarios",
]
