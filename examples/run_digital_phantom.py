"""Run software verification using the 1 mV / 60 bpm digital calibration phantom.

Usage:
    python examples/run_digital_phantom.py
"""

from __future__ import annotations

import json

from ecg_integrity.digital_phantom import evaluate_digital_phantom, generate_1mv_60bpm_phantom


def main() -> None:
    # 1. Nominal 1.0 mV, 60 bpm, 360 Hz, 10 s, 100 ms pulse width verification
    phantom_nominal = generate_1mv_60bpm_phantom(
        duration_s=10.0,
        fs=360.0,
        amplitude_mv=1.0,
        bpm=60.0,
        pulse_width_ms=100.0,
        gain_factor=1.0,
    )
    nominal_report = evaluate_digital_phantom(
        phantom_nominal,
        fs=360.0,
        expected_amplitude_mv=1.0,
        expected_bpm=60.0,
    )
    print("=== Nominal Digital Calibration Phantom (1.0 mV, 60 bpm, 360 Hz) ===")
    print(json.dumps(nominal_report, indent=2))

    # 2. Sensitivity verification under an artificial +5% scale deviation (1.05 mV)
    phantom_gain_5pct = generate_1mv_60bpm_phantom(
        duration_s=10.0,
        fs=360.0,
        amplitude_mv=1.0,
        bpm=60.0,
        pulse_width_ms=100.0,
        gain_factor=1.05,
    )
    gain_report = evaluate_digital_phantom(
        phantom_gain_5pct,
        fs=360.0,
        expected_amplitude_mv=1.0,
        expected_bpm=60.0,
    )
    print("\n=== +5% Gain Perturbation Check (1.05 mV vs 1.0 mV reference) ===")
    print(json.dumps(gain_report, indent=2))


if __name__ == "__main__":
    main()
