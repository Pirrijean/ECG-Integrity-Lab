"""PhysioNet MIT-BIH Noise Stress Test Database (nstdb) loader and SNR mixer.

Implements Experiment 2 (Real Physiological Acquisition Noise):
Loads physical noise recordings collected from active volunteers:
- ``bw``: baseline wander
- ``ma``: muscle (EMG) artifact
- ``em``: electrode motion artifact

Mixes noise ``n[n]`` onto clean calibrated ECG ``x[n]`` at controlled target SNRs
(e.g., 24, 18, 12, 6, 0, -6 dB) while preserving exact beat annotation alignment:
    x_noisy[n] = x[n] + alpha * n_0[n]
    alpha = sqrt( P_x / (P_n * 10^(SNR_dB / 10)) )
where P_x = Var(x) and P_n = Var(n_0) are computed over the aligned segment.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import wfdb

from ecg_integrity.signal_processing import validate_signal_array

NSTDB_NOISE_TYPES: tuple[str, ...] = ("bw", "ma", "em")
NSTDB_DEFAULT_SNRS_DB: tuple[float, ...] = (24.0, 18.0, 12.0, 6.0, 0.0, -6.0)

NSTDB_LABELS: dict[str, str] = {
    "bw": "Real Baseline Wander (nstdb/bw)",
    "ma": "Real Muscle Artifact (nstdb/ma)",
    "em": "Real Electrode Motion (nstdb/em)",
}


def load_nstdb_noise(
    noise_type: str,
    n_samples: int,
    channel: int = 0,
    expected_fs: float = 360.0,
    pn_dir: str = "nstdb",
    cache_dir: str | Path | None = "data/raw/nstdb",
    offset_samples: int = 0,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Load and validate a physical noise segment from PhysioNet's nstdb."""
    ntype = noise_type.strip().lower()
    if ntype not in NSTDB_NOISE_TYPES:
        raise ValueError(f"noise_type must be one of {NSTDB_NOISE_TYPES}; got {noise_type!r}.")
    if n_samples < 2:
        raise ValueError("n_samples must be >= 2.")

    local_path: str | None = None
    if cache_dir is not None:
        cpath = Path(cache_dir)
        cpath.mkdir(parents=True, exist_ok=True)
        req = [f"{ntype}.hea", f"{ntype}.dat"]
        if not all((cpath / f).exists() for f in req):
            wfdb.dl_files(pn_dir, str(cpath), req)
        local_path = str(cpath / ntype)

    sampto = offset_samples + n_samples
    if local_path is not None:
        rec = wfdb.rdrecord(local_path, sampfrom=offset_samples, sampto=sampto, physical=True)
    else:
        rec = wfdb.rdrecord(ntype, sampfrom=offset_samples, sampto=sampto, physical=True, pn_dir=pn_dir)

    fs = float(rec.fs)
    if abs(fs - expected_fs) > 1e-6:
        raise ValueError(f"NSTDB {ntype} fs={fs} Hz does not match expected ECG fs={expected_fs} Hz.")
    if rec.p_signal is None or channel < 0 or channel >= rec.p_signal.shape[1]:
        raise ValueError(f"Invalid channel {channel} on NSTDB record {ntype}.")

    unit = str(rec.units[channel]).strip()
    if unit.lower() != "mv":
        raise ValueError(f"Expected mV physical unit on NSTDB {ntype}, got {unit!r}.")

    raw_noise = np.ascontiguousarray(rec.p_signal[:, channel], dtype=np.float64)
    if raw_noise.size < n_samples:
        # Tile seamlessly only if requested duration exceeds the 30-min noise tape
        reps = int(np.ceil(n_samples / raw_noise.size))
        raw_noise = np.tile(raw_noise, reps)[:n_samples]

    if not np.all(np.isfinite(raw_noise)):
        raise ValueError(f"NSTDB {ntype} contains non-finite samples.")

    meta = {
        "noise_type": ntype,
        "label": NSTDB_LABELS[ntype],
        "fs_hz": fs,
        "channel": channel,
        "unit": unit,
        "n_samples": int(raw_noise.size),
        "raw_variance_mv2": round(float(np.var(raw_noise)), 8),
    }
    return raw_noise, meta


def mix_ecg_with_real_noise(
    clean_ecg: np.ndarray,
    noise_segment: np.ndarray,
    snr_db: float,
    fs: float = 360.0,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Combine clean ECG x[n] with zero-mean physical noise n[n] at target SNR (dB).

    Uses segment AC power definitions:
        P_x = Var(x),   P_n = Var(n - mean(n))
        alpha = sqrt( P_x / (P_n * 10^(SNR_dB / 10)) )
    """
    x = validate_signal_array(clean_ecg, fs=fs, max_abs_mv=50.0)
    n_arr = np.asarray(noise_segment, dtype=np.float64)
    if n_arr.ndim != 1 or n_arr.size != x.size:
        raise ValueError(f"noise_segment length ({n_arr.size}) must equal clean_ecg length ({x.size}).")
    if not np.all(np.isfinite(n_arr)):
        raise ValueError("noise_segment contains non-finite values.")
    if not np.isfinite(snr_db):
        raise ValueError(f"snr_db must be finite; got {snr_db!r}.")

    n_zero = n_arr - float(np.mean(n_arr))
    p_x = float(np.var(x))
    p_n = float(np.var(n_zero))
    if p_x <= 1e-15 or p_n <= 1e-15:
        raise ValueError("Signal power P_x and noise power P_n must both be strictly positive.")

    alpha = float(np.sqrt(p_x / (p_n * (10.0 ** (float(snr_db) / 10.0)))))
    scaled_noise = alpha * n_zero
    mixed = x + scaled_noise
    empirical_snr = float(10.0 * np.log10(p_x / float(np.var(scaled_noise))))

    mix_info = {
        "target_snr_db": float(snr_db),
        "empirical_snr_db": round(empirical_snr, 4),
        "signal_variance_mv2": round(p_x, 8),
        "noise_variance_mv2": round(p_n, 8),
        "alpha_scale": round(alpha, 6),
    }
    return mixed, mix_info
