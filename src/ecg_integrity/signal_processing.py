"""ECG signal conditioning (Butterworth SOS zero-phase filter) and controlled perturbations.

Scope warning:
The 5–20 Hz bandpass filter implemented here is designed solely to emphasize QRS
energy for beat localization. It attenuates P/T waves and can distort ST-segment
morphology; do NOT use its output for diagnostic waveform or ST-segment interpretation.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.signal import butter, sosfiltfilt


def validate_signal_array(
    signal: np.ndarray,
    fs: float,
    max_abs_mv: float = 50.0,
) -> np.ndarray:
    """Validate 1D ECG signal array, finiteness, amplitude bounds, and sampling rate."""
    if not np.isfinite(fs) or fs <= 0.0:
        raise ValueError(f"Sampling frequency 'fs' must be finite and > 0 Hz; got {fs!r}.")

    arr = np.asarray(signal, dtype=np.float64)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError("ECG signal must be a non-empty 1D array.")
    if not np.all(np.isfinite(arr)):
        raise ValueError("ECG signal contains non-finite values (NaN or Inf).")
    if np.max(np.abs(arr)) > max_abs_mv:
        raise ValueError(
            f"Signal amplitude ({float(np.max(np.abs(arr))):.3f} mV) exceeds allowed bound ({max_abs_mv:.3f} mV)."
        )
    return arr


def bandpass_filter(
    signal: np.ndarray,
    fs: float,
    lowcut_hz: float = 5.0,
    highcut_hz: float = 20.0,
    order: int = 4,
) -> np.ndarray:
    """Apply a zero-phase Butterworth bandpass filter using second-order sections (SOS).

    Uses ``scipy.signal.butter(..., output='sos')`` and ``scipy.signal.sosfiltfilt``
    for numerical stability and zero phase distortion of QRS peak timing.
    """
    x = validate_signal_array(signal, fs)
    nyquist = 0.5 * fs
    if not (0.0 < lowcut_hz < highcut_hz < nyquist):
        raise ValueError(
            f"Cutoff frequencies must satisfy 0 < lowcut ({lowcut_hz}) < highcut ({highcut_hz}) < Nyquist ({nyquist})."
        )
    if not isinstance(order, (int, np.integer)) or order < 1:
        raise ValueError(f"Filter order must be a positive integer; got {order!r}.")

    sos = butter(int(order), [lowcut_hz, highcut_hz], btype="bandpass", fs=fs, output="sos")
    return np.asarray(sosfiltfilt(sos, x), dtype=np.float64)


def add_baseline_wander(
    signal: np.ndarray,
    fs: float,
    amplitude_mv: float = 0.20,
    freq_hz: float = 0.25,
    phase_rad: float = 0.0,
) -> np.ndarray:
    """Superimpose sinusoidal respiratory/motion baseline wander onto an ECG signal."""
    x = validate_signal_array(signal, fs)
    if not np.isfinite(amplitude_mv) or amplitude_mv < 0.0:
        raise ValueError("amplitude_mv must be finite and non-negative.")
    if not np.isfinite(freq_hz) or freq_hz <= 0.0:
        raise ValueError("freq_hz must be finite and positive.")

    t = np.arange(x.size, dtype=np.float64) / fs
    wander = amplitude_mv * np.sin(2.0 * np.pi * freq_hz * t + phase_rad)
    return x + wander


def add_mains_interference(
    signal: np.ndarray,
    fs: float,
    amplitude_mv: float = 0.05,
    freq_hz: float = 60.0,
    phase_rad: float = 0.0,
) -> np.ndarray:
    """Superimpose sinusoidal power-line interference (e.g., 60 Hz) onto an ECG signal."""
    x = validate_signal_array(signal, fs)
    if not np.isfinite(amplitude_mv) or amplitude_mv < 0.0:
        raise ValueError("amplitude_mv must be finite and non-negative.")
    if not np.isfinite(freq_hz) or freq_hz <= 0.0 or freq_hz >= 0.5 * fs:
        raise ValueError(f"Mains frequency ({freq_hz} Hz) must lie in (0, Nyquist={0.5 * fs} Hz).")

    t = np.arange(x.size, dtype=np.float64) / fs
    mains = amplitude_mv * np.sin(2.0 * np.pi * freq_hz * t + phase_rad)
    return x + mains


def add_white_noise(
    signal: np.ndarray,
    snr_db: float,
    seed: int = 42,
    fs: float = 360.0,
) -> np.ndarray:
    """Add reproducible zero-mean Gaussian broadband noise at a specified SNR (dB).

    Signal power is defined on the AC (zero-mean) waveform:
        P_signal = mean((x - mean(x))^2)
        P_noise  = P_signal / 10^(SNR_dB / 10)
    so that DC baseline offsets do not artificially inflate the estimated signal power.
    """
    x = validate_signal_array(signal, fs, max_abs_mv=100.0)
    if not np.isfinite(snr_db):
        raise ValueError(f"snr_db must be finite; got {snr_db!r}.")

    sig_power = float(np.var(x))
    if sig_power <= 1e-15:
        sig_power = float(np.mean(x**2))
    if sig_power <= 1e-15:
        raise ValueError("Cannot compute SNR-scaled noise on an identically zero signal.")

    noise_power = sig_power / (10.0 ** (float(snr_db) / 10.0))
    noise_std = float(np.sqrt(noise_power))

    rng = np.random.default_rng(seed)
    noise = rng.normal(loc=0.0, scale=noise_std, size=x.shape)
    return x + noise


def build_protocol_scenarios(
    clean_signal: np.ndarray,
    fs: float,
    seed: int = 42,
) -> list[dict[str, Any]]:
    """Construct the 6 experimental interference scenarios defined in Step 15."""
    x = validate_signal_array(clean_signal, fs)
    return [
        {
            "key": "baseline",
            "scenario": "Baseline",
            "condition": "Original signal",
            "signal": x.copy(),
            "parameters": {"perturbation": "none"},
        },
        {
            "key": "baseline_wander",
            "scenario": "Baseline wander",
            "condition": "0.20 mV, 0.25 Hz",
            "signal": add_baseline_wander(x, fs, amplitude_mv=0.20, freq_hz=0.25),
            "parameters": {"type": "sinusoidal_wander", "amplitude_mv": 0.20, "freq_hz": 0.25},
        },
        {
            "key": "mains_interference",
            "scenario": "Mains interference",
            "condition": "0.05 mV, 60 Hz",
            "signal": add_mains_interference(x, fs, amplitude_mv=0.05, freq_hz=60.0),
            "parameters": {"type": "powerline_60hz", "amplitude_mv": 0.05, "freq_hz": 60.0},
        },
        {
            "key": "white_noise_20dB",
            "scenario": "White noise (20 dB)",
            "condition": "AWGN 20 dB SNR",
            "signal": add_white_noise(x, snr_db=20.0, seed=seed, fs=fs),
            "parameters": {"type": "awgn", "snr_db": 20.0, "seed": seed},
        },
        {
            "key": "white_noise_10dB",
            "scenario": "White noise (10 dB)",
            "condition": "AWGN 10 dB SNR",
            "signal": add_white_noise(x, snr_db=10.0, seed=seed, fs=fs),
            "parameters": {"type": "awgn", "snr_db": 10.0, "seed": seed},
        },
        {
            "key": "white_noise_5dB",
            "scenario": "White noise (5 dB)",
            "condition": "AWGN 5 dB SNR",
            "signal": add_white_noise(x, snr_db=5.0, seed=seed, fs=fs),
            "parameters": {"type": "awgn", "snr_db": 5.0, "seed": seed},
        },
    ]
