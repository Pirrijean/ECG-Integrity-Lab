"""Unit and QA verification suite for ECG-Integrity-Lab."""

from __future__ import annotations

import numpy as np
import pytest

from ecg_integrity.dataset import PHYSIONET_BEAT_SYMBOLS, validate_ecg_metadata_and_signal
from ecg_integrity.detector import detect_qrs
from ecg_integrity.digital_phantom import evaluate_digital_phantom, generate_1mv_60bpm_phantom
from ecg_integrity.evaluation import match_beats_one_to_one
from ecg_integrity.signal_processing import add_white_noise, bandpass_filter


def test_one_to_one_beat_matching_prevents_double_counting() -> None:
    """Two detections near a single reference beat must yield 1 TP and 1 FP, never 2 TPs."""
    fs = 360.0
    ref_peaks = np.array([360, 720, 1080], dtype=np.int64)
    # Near 360: both 362 (+2 samples) and 380 (+20 samples) are within 150 ms (54 samples).
    # 720 is missed (FN). 1080 is matched by 1078 (-2 samples).
    det_peaks = np.array([362, 380, 1078], dtype=np.int64)

    res = match_beats_one_to_one(ref_peaks, det_peaks, fs=fs, tolerance_ms=150.0)
    assert res.tp == 2
    assert res.fp == 1
    assert res.fn == 1
    assert list(res.matched_ref_samples) == [360, 1080]
    assert list(res.matched_det_samples) == [362, 1078]
    assert list(res.unmatched_det_samples) == [380]
    assert list(res.unmatched_ref_samples) == [720]


def test_evaluation_metrics_sensitivity_ppv_f1_and_mae() -> None:
    """Verify exact analytical formulas for Sensitivity, PPV, F1-score, and MAE_t (ms)."""
    fs = 360.0
    ref_peaks = np.array([100, 500, 900, 1300], dtype=np.int64)
    # 3 TPs with sample offsets +1, -1, 0; 1 FN (1300); 1 FP (1600)
    det_peaks = np.array([101, 499, 900, 1600], dtype=np.int64)

    res = match_beats_one_to_one(ref_peaks, det_peaks, fs=fs, tolerance_ms=150.0)
    assert res.tp == 3
    assert res.fp == 1
    assert res.fn == 1
    assert res.sensitivity == pytest.approx(0.75)
    assert res.ppv == pytest.approx(0.75)
    assert res.f1_score == pytest.approx(0.75)
    expected_mae_ms = ((1 + 1 + 0) / 3.0) * (1000.0 / 360.0)
    assert res.mae_ms == pytest.approx(expected_mae_ms, rel=1e-6)


def test_butterworth_bandpass_attenuates_wander_and_mains_preserving_qrs_band() -> None:
    """Verify 5–20 Hz Butterworth SOS filter suppresses 0.25 Hz and 60 Hz while passing 12 Hz."""
    fs = 360.0
    t = np.arange(int(20 * fs), dtype=np.float64) / fs
    wander = 0.50 * np.sin(2.0 * np.pi * 0.25 * t)
    qrs_band = 1.00 * np.sin(2.0 * np.pi * 12.0 * t)
    mains = 0.50 * np.sin(2.0 * np.pi * 60.0 * t)

    composite = wander + qrs_band + mains
    filtered = bandpass_filter(composite, fs=fs, lowcut_hz=5.0, highcut_hz=20.0, order=4)

    # Trim 2 seconds at edges to ignore transient edge effects
    mid = slice(int(2 * fs), int(18 * fs))
    residual = filtered[mid] - qrs_band[mid]
    assert np.std(residual) < 0.02
    assert np.corrcoef(filtered[mid], qrs_band[mid])[0, 1] > 0.999


def test_white_noise_reproducibility_and_snr_scaling() -> None:
    """Same random seed must produce identical noisy waveforms and accurate empirical SNR."""
    fs = 360.0
    t = np.arange(int(30 * fs), dtype=np.float64) / fs
    clean = 0.8 * np.sin(2.0 * np.pi * 10.0 * t) + 0.35  # Includes DC offset

    noisy_a = add_white_noise(clean, snr_db=10.0, seed=42, fs=fs)
    noisy_b = add_white_noise(clean, snr_db=10.0, seed=42, fs=fs)
    noisy_c = add_white_noise(clean, snr_db=10.0, seed=99, fs=fs)

    assert np.array_equal(noisy_a, noisy_b)
    assert not np.array_equal(noisy_a, noisy_c)

    injected_noise = noisy_a - clean
    empirical_snr_db = 10.0 * np.log10(np.var(clean) / np.var(injected_noise))
    assert empirical_snr_db == pytest.approx(10.0, abs=0.25)


def test_digital_calibration_phantom_nominal_and_gain_error() -> None:
    """Verify 0.0% error on 1 mV / 60 bpm phantom and exact +5.0% error on 1.05x gain."""
    phantom = generate_1mv_60bpm_phantom(duration_s=10.0, fs=360.0, amplitude_mv=1.0, bpm=60.0)
    report = evaluate_digital_phantom(phantom, fs=360.0, expected_amplitude_mv=1.0, expected_bpm=60.0)

    assert report["pulse_count"] == 10
    assert report["amplitude_error_percent"] == pytest.approx(0.0, abs=1e-6)
    assert report["rate_error_percent"] == pytest.approx(0.0, abs=1e-6)

    phantom_5pct = generate_1mv_60bpm_phantom(
        duration_s=10.0, fs=360.0, amplitude_mv=1.0, bpm=60.0, gain_factor=1.05
    )
    report_5pct = evaluate_digital_phantom(
        phantom_5pct, fs=360.0, expected_amplitude_mv=1.0, expected_bpm=60.0
    )
    assert report_5pct["observed_amplitude_mv"] == pytest.approx(1.05, abs=1e-6)
    assert report_5pct["amplitude_error_percent"] == pytest.approx(5.0, abs=1e-6)


def test_qa_validation_rejects_invalid_metadata_and_corrupted_signals() -> None:
    """Section 11.A: verify rejection of invalid fs, NaN/Inf, nonexistent channel, and out-of-range values."""
    good_sig = np.zeros((3600, 2), dtype=np.float64)
    good_sig[180::360, 0] = 1.0
    good_refs = np.arange(180, 3600, 360, dtype=np.int64)

    # 1. Invalid sampling frequency
    with pytest.raises(ValueError, match="sampling frequency"):
        validate_ecg_metadata_and_signal(good_sig, 0, fs=0.0, units_list=["mV", "mV"], sig_names=["MLII", "V5"], r_peaks_ref=good_refs)

    # 2. Nonexistent channel index
    with pytest.raises(ValueError, match="channel"):
        validate_ecg_metadata_and_signal(good_sig, 5, fs=360.0, units_list=["mV", "mV"], sig_names=["MLII", "V5"], r_peaks_ref=good_refs)

    # 3. Signal containing NaN
    nan_sig = good_sig.copy()
    nan_sig[50, 0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        validate_ecg_metadata_and_signal(nan_sig, 0, fs=360.0, units_list=["mV", "mV"], sig_names=["MLII", "V5"], r_peaks_ref=good_refs)

    # 4. Out-of-range uncalibrated ADC-like values (e.g. 1024 counts instead of mV)
    raw_adc_sig = np.full((3600, 2), 1024.0, dtype=np.float64)
    with pytest.raises(ValueError, match="range"):
        validate_ecg_metadata_and_signal(raw_adc_sig, 0, fs=360.0, units_list=["mV", "mV"], sig_names=["MLII", "V5"], r_peaks_ref=good_refs)

    # 5. Beat symbol filter excludes rhythm/artifact markers
    assert "N" in PHYSIONET_BEAT_SYMBOLS
    assert "V" in PHYSIONET_BEAT_SYMBOLS
    assert "+" not in PHYSIONET_BEAT_SYMBOLS
    assert "~" not in PHYSIONET_BEAT_SYMBOLS


def test_detector_recovers_synthetic_qrs_train() -> None:
    """Verify Pan–Tompkins detector locates synthetic Gaussian QRS complexes accurately."""
    fs = 360.0
    n = int(12 * fs)
    t = np.arange(n, dtype=np.float64) / fs
    sig = np.zeros(n, dtype=np.float64)
    true_peaks = np.arange(int(0.8 * fs), n - int(0.5 * fs), int(0.85 * fs), dtype=np.int64)
    for p in true_peaks:
        tp = p / fs
        sig += 1.2 * np.exp(-0.5 * ((t - tp) / 0.012) ** 2)

    detected = detect_qrs(sig, fs=fs)
    eval_res = match_beats_one_to_one(true_peaks, detected, fs=fs, tolerance_ms=50.0)
    assert eval_res.sensitivity == pytest.approx(1.0)
    assert eval_res.ppv == pytest.approx(1.0)
    assert eval_res.mae_ms < 5.0
