"""Public package interface for ECG-Integrity-Lab."""

from __future__ import annotations

__version__ = "0.1.0"

from ecg_integrity.dataset import (
    PHYSIONET_BEAT_SYMBOLS,
    ECGRecordData,
    load_mitbih_record,
    validate_ecg_metadata_and_signal,
)
from ecg_integrity.detector import QRSDetectionOutput, QRSDetectorConfig, detect_qrs, detect_qrs_pan_tompkins
from ecg_integrity.digital_phantom import evaluate_digital_phantom, generate_1mv_60bpm_phantom
from ecg_integrity.evaluation import BeatEvaluationResult, match_beats_one_to_one
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
    "QRSDetectionOutput",
    "detect_qrs",
    "detect_qrs_pan_tompkins",
    "generate_1mv_60bpm_phantom",
    "evaluate_digital_phantom",
    "BeatEvaluationResult",
    "match_beats_one_to_one",
    "bandpass_filter",
    "add_baseline_wander",
    "add_mains_interference",
    "add_white_noise",
    "build_protocol_scenarios",
]
