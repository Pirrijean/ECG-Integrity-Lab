"""Functional and mathematical software verification via digital reference signals.

Important engineering note:
This module generates a deterministic digital pulse train to verify software
amplitude scaling and temporal interval arithmetic (functional software test).
It does NOT verify an electrocardiograph's physical acquisition chain, is NOT
traceable to SI units, and does NOT constitute metrological hardware calibration.
"""

from __future__ import annotations

from typing import Any

import numpy as np


def generate_1mv_60bpm_phantom(
    duration_s: float = 10.0,
    fs: float = 360.0,
    amplitude_mv: float = 1.0,
    bpm: float = 60.0,
    pulse_width_ms: float = 100.0,
    gain_factor: float = 1.0,
) -> np.ndarray:
    """Generate a synthetic rectangular pulse train with known amplitude and rate.

    Parameters
    ----------
    duration_s : float, default=10.0
        Signal duration in seconds.
    fs : float, default=360.0
        Sampling frequency in Hz (dt = 1 / fs).
    amplitude_mv : float, default=1.0
        Nominal pulse height in millivolts.
    bpm : float, default=60.0
        Nominal pulse repetition rate in beats per minute.
    pulse_width_ms : float, default=100.0
        Width of each rectangular pulse in milliseconds.
    gain_factor : float, default=1.0
        Multiplicative scale factor (e.g., 1.05 to inject +5% gain error).

    Returns
    -------
    np.ndarray
        1D float64 array representing the synthetic phantom signal in mV.
    """
    for name, val in (
        ("duration_s", duration_s),
        ("fs", fs),
        ("amplitude_mv", amplitude_mv),
        ("bpm", bpm),
        ("pulse_width_ms", pulse_width_ms),
        ("gain_factor", gain_factor),
    ):
        if not np.isfinite(val) or val <= 0.0:
            raise ValueError(f"Parameter '{name}' must be finite and strictly positive; got {val!r}.")

    period_s = 60.0 / bpm
    pulse_width_s = pulse_width_ms / 1000.0
    if pulse_width_s >= period_s:
        raise ValueError(
            f"Pulse width ({pulse_width_s:.4f} s) must be strictly smaller than beat period ({period_s:.4f} s)."
        )

    n_samples = int(round(duration_s * fs))
    if n_samples < 2:
        raise ValueError("Signal duration and sampling frequency yield fewer than 2 samples.")

    signal = np.zeros(n_samples, dtype=np.float64)
    period_samples = period_s * fs
    width_samples = max(1, int(round(pulse_width_s * fs)))
    half_width = width_samples // 2

    # Place pulses centered at (k + 0.5) * period_samples so all pulses lie cleanly inside [0, duration_s)
    k = 0
    target_amplitude = float(amplitude_mv * gain_factor)
    while True:
        center = int(round((k + 0.5) * period_samples))
        start = center - half_width
        end = start + width_samples
        if end > n_samples:
            break
        if start >= 0:
            signal[start:end] = target_amplitude
        k += 1

    return signal


def evaluate_digital_phantom(
    signal: np.ndarray,
    fs: float = 360.0,
    expected_amplitude_mv: float = 1.0,
    expected_bpm: float = 60.0,
) -> dict[str, Any]:
    """Verify amplitude and rate accuracy of a digital calibration pulse train.

    Computes percentage scale error:
        E_A = ((A_measured - A_reference) / A_reference) * 100
    and percentage rate error:
        E_R = ((BPM_measured - BPM_reference) / BPM_reference) * 100
    """
    arr = np.asarray(signal, dtype=np.float64)
    if arr.ndim != 1 or arr.size < 2:
        raise ValueError("Phantom signal must be a non-empty 1D array with at least 2 samples.")
    if not np.all(np.isfinite(arr)):
        raise ValueError("Phantom signal contains non-finite values (NaN or Inf).")
    if not np.isfinite(fs) or fs <= 0.0:
        raise ValueError(f"Sampling frequency 'fs' must be finite and positive; got {fs!r}.")
    if not np.isfinite(expected_amplitude_mv) or expected_amplitude_mv <= 0.0:
        raise ValueError("expected_amplitude_mv must be finite and positive.")
    if not np.isfinite(expected_bpm) or expected_bpm <= 0.0:
        raise ValueError("expected_bpm must be finite and positive.")

    baseline_mv = float(np.min(arr))
    peak_mv = float(np.max(arr))
    dynamic_range = peak_mv - baseline_mv
    if dynamic_range <= 0.0:
        raise ValueError("Phantom signal is flat; no pulses detected.")

    threshold = baseline_mv + 0.5 * dynamic_range
    active = arr > threshold

    # Identify rising and falling transitions with boolean padding
    padded = np.concatenate(([False], active, [False]))
    transitions = np.diff(padded.astype(np.int8))
    starts = np.flatnonzero(transitions == 1)
    ends = np.flatnonzero(transitions == -1)

    pulse_count = int(starts.size)
    if pulse_count < 2:
        raise ValueError(f"At least 2 pulses are required to compute rate; detected {pulse_count}.")

    pulse_heights = [float(np.max(arr[s:e]) - baseline_mv) for s, e in zip(starts, ends)]
    pulse_centers = np.array([0.5 * (s + e - 1) for s, e in zip(starts, ends)], dtype=np.float64)

    observed_amplitude_mv = float(np.mean(pulse_heights))
    amplitude_error_percent = ((observed_amplitude_mv - expected_amplitude_mv) / expected_amplitude_mv) * 100.0

    rr_intervals_s = np.diff(pulse_centers) / fs
    mean_rr_s = float(np.mean(rr_intervals_s))
    observed_bpm = 60.0 / mean_rr_s
    rate_error_percent = ((observed_bpm - expected_bpm) / expected_bpm) * 100.0

    return {
        "pulse_count": pulse_count,
        "expected_amplitude_mv": round(float(expected_amplitude_mv), 6),
        "observed_amplitude_mv": round(observed_amplitude_mv, 6),
        "amplitude_error_percent": round(float(amplitude_error_percent), 6),
        "expected_bpm": round(float(expected_bpm), 6),
        "observed_bpm": round(float(observed_bpm), 6),
        "rate_error_percent": round(float(rate_error_percent), 6),
    }
