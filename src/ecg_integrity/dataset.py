"""PhysioNet / MIT-BIH Arrhythmia Database loader and metadata QA validator."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import wfdb

# Standard PhysioNet beat annotation symbols (true QRS complexes).
# Excludes rhythm markers ('+'), signal-quality markers ('~'), isolated artifacts ('|'),
# non-conducted P-waves ('x'), waveform boundary markers, and text notes.
PHYSIONET_BEAT_SYMBOLS: frozenset[str] = frozenset(
    {
        "N",  # Normal beat
        "L",  # Left bundle branch block beat
        "R",  # Right bundle branch block beat
        "B",  # Bundle branch block beat (unspecified)
        "A",  # Atrial premature beat
        "a",  # Aberrated atrial premature beat
        "J",  # Nodal (junctional) premature beat
        "S",  # Supraventricular premature or ectopic beat
        "V",  # Premature ventricular contraction
        "r",  # R-on-T premature ventricular contraction
        "F",  # Fusion of ventricular and normal beat
        "e",  # Atrial escape beat
        "j",  # Nodal (junctional) escape beat
        "n",  # Supraventricular escape beat
        "E",  # Ventricular escape beat
        "/",  # Paced beat
        "f",  # Fusion of paced and normal beat
        "Q",  # Unclassifiable beat
        "?",  # Beat not classified during learning
    }
)


@dataclass(frozen=True)
class ECGRecordData:
    """Validated ECG signal segment, expert beat annotations, and header metadata."""

    record_id: str
    channel: int
    channel_name: str
    fs: float
    units: str
    signal: np.ndarray
    r_peaks_ref: np.ndarray
    beat_symbols: list[str]
    all_annotation_count: int
    non_beat_annotation_count: int
    metadata: dict[str, Any]


def validate_ecg_metadata_and_signal(
    signal_matrix: np.ndarray,
    channel: int,
    fs: float,
    units_list: list[str],
    sig_names: list[str],
    r_peaks_ref: np.ndarray | None = None,
    expected_fs: float | None = 360.0,
    expected_unit: str = "mV",
    max_abs_mv: float = 20.0,
    require_annotations: bool = True,
    adc_gain: float | None = None,
    adc_baseline: int | None = None,
    adc_res: int | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Perform Step 9 QA/metrology verification on header metadata and physical samples.

    Raises ValueError on invalid sampling rate, nonexistent channel, wrong physical
    units, non-finite values (NaN/Inf), out-of-range amplitudes, or missing annotations.
    """
    if not np.isfinite(fs) or fs <= 0.0:
        raise ValueError(f"Invalid sampling frequency fs={fs!r}; must be finite and > 0 Hz.")

    if expected_fs is not None and abs(float(fs) - float(expected_fs)) > 1e-6:
        raise ValueError(
            f"Sampling frequency mismatch: expected {expected_fs} Hz, observed {fs} Hz."
        )

    arr2d = np.asarray(signal_matrix, dtype=np.float64)
    if arr2d.ndim == 1:
        arr2d = arr2d[:, np.newaxis]
    if arr2d.ndim != 2 or arr2d.shape[0] == 0:
        raise ValueError("Signal array must be non-empty 1D or 2D.")

    n_samples, n_channels = arr2d.shape
    if not isinstance(channel, (int, np.integer)) or channel < 0 or channel >= n_channels:
        raise ValueError(
            f"Requested channel {channel} does not exist (record has {n_channels} channel(s): 0..{n_channels - 1})."
        )

    if len(units_list) <= channel:
        raise ValueError(f"Missing physical unit entry for channel {channel}.")
    unit = str(units_list[channel]).strip()
    if unit.lower() != expected_unit.strip().lower():
        raise ValueError(
            f"Physical unit mismatch on channel {channel}: expected {expected_unit!r}, got {unit!r}."
        )

    channel_signal = np.ascontiguousarray(arr2d[:, channel], dtype=np.float64)
    if not np.all(np.isfinite(channel_signal)):
        n_bad = int(np.sum(~np.isfinite(channel_signal)))
        raise ValueError(f"Signal contains {n_bad} non-finite sample(s) (NaN or Inf).")

    max_observed_abs = float(np.max(np.abs(channel_signal)))
    if max_observed_abs > max_abs_mv:
        raise ValueError(
            f"Signal amplitude out of physiological/calibrated range: max |x| = {max_observed_abs:.3f} {unit} "
            f"exceeds limit of {max_abs_mv:.3f} {expected_unit}. Verify ADC calibration scaling."
        )

    if require_annotations:
        if r_peaks_ref is None or len(r_peaks_ref) == 0:
            raise ValueError("No valid reference beat annotations available in the selected window.")

    channel_name = str(sig_names[channel]) if len(sig_names) > channel else f"ch{channel}"
    duration_s = float(n_samples) / float(fs)
    dt_s = 1.0 / float(fs)

    qa_metadata: dict[str, Any] = {
        "n_channels": int(n_channels),
        "selected_channel": int(channel),
        "channel_name": channel_name,
        "fs_hz": float(fs),
        "expected_fs_hz": float(expected_fs) if expected_fs is not None else None,
        "dt_seconds": round(dt_s, 8),
        "dt_ms": round(dt_s * 1000.0, 6),
        "physical_unit": unit,
        "n_samples": int(n_samples),
        "duration_seconds": round(duration_s, 4),
        "non_finite_count": 0,
        "min_amplitude_mv": round(float(np.min(channel_signal)), 6),
        "max_amplitude_mv": round(float(np.max(channel_signal)), 6),
        "mean_amplitude_mv": round(float(np.mean(channel_signal)), 6),
        "std_amplitude_mv": round(float(np.std(channel_signal)), 6),
        "adc_gain_counts_per_mv": float(adc_gain) if adc_gain is not None else None,
        "adc_baseline_counts": int(adc_baseline) if adc_baseline is not None else None,
        "adc_resolution_bits": int(adc_res) if adc_res is not None else None,
        "reference_beat_count": int(len(r_peaks_ref)) if r_peaks_ref is not None else 0,
        "qa_checks_passed": True,
    }
    return channel_signal, qa_metadata


def load_mitbih_record(
    record_id: str | int,
    channel: int = 0,
    seconds: float | None = 300.0,
    pn_dir: str = "mitdb",
    cache_dir: str | Path | None = "data/raw/mitdb",
    expected_fs: float | None = 360.0,
) -> ECGRecordData:
    """Load an ECG record and expert '.atr' beat annotations from MIT-BIH via WFDB.

    Caches header (.hea), signal (.dat), and annotation (.atr) files locally on
    first access so subsequent runs are deterministic and fast.
    """
    rec_str = str(record_id).strip()
    if not rec_str:
        raise ValueError("record_id must be a non-empty string or integer.")
    if seconds is not None and (not np.isfinite(seconds) or seconds <= 0.0):
        raise ValueError(f"seconds must be finite and > 0 when provided; got {seconds!r}.")

    local_record_path: str | None = None
    if cache_dir is not None:
        cache_path = Path(cache_dir)
        cache_path.mkdir(parents=True, exist_ok=True)
        required_files = [f"{rec_str}.hea", f"{rec_str}.dat", f"{rec_str}.atr"]
        if not all((cache_path / fname).exists() for fname in required_files):
            wfdb.dl_files(pn_dir, str(cache_path), required_files)
        local_record_path = str(cache_path / rec_str)

    if local_record_path is not None:
        header = wfdb.rdheader(local_record_path)
    else:
        header = wfdb.rdheader(rec_str, pn_dir=pn_dir)

    fs_header = float(header.fs)
    total_samples = int(header.sig_len)
    sampto = total_samples
    if seconds is not None:
        requested_samples = int(round(seconds * fs_header))
        sampto = min(total_samples, requested_samples)

    if local_record_path is not None:
        record = wfdb.rdrecord(local_record_path, sampfrom=0, sampto=sampto, physical=True)
        ann = wfdb.rdann(local_record_path, "atr", sampfrom=0, sampto=sampto)
    else:
        record = wfdb.rdrecord(rec_str, sampfrom=0, sampto=sampto, physical=True, pn_dir=pn_dir)
        ann = wfdb.rdann(rec_str, "atr", sampfrom=0, sampto=sampto, pn_dir=pn_dir)

    if record.p_signal is None:
        raise ValueError(f"Record {rec_str} did not return calibrated physical signals (p_signal is None).")

    # Filter PhysioNet annotations strictly to true QRS beat symbols within [0, sampto)
    beat_samples: list[int] = []
    beat_syms: list[str] = []
    non_beat_count = 0
    for samp, sym in zip(ann.sample, ann.symbol):
        samp_int = int(samp)
        if 0 <= samp_int < record.p_signal.shape[0]:
            if sym in PHYSIONET_BEAT_SYMBOLS:
                beat_samples.append(samp_int)
                beat_syms.append(str(sym))
            else:
                non_beat_count += 1

    r_peaks_ref = np.asarray(beat_samples, dtype=np.int64)

    adc_gain = None
    if getattr(record, "adc_gain", None) is not None and len(record.adc_gain) > channel >= 0:
        adc_gain = float(record.adc_gain[channel])
    adc_baseline = None
    if getattr(record, "baseline", None) is not None and len(record.baseline) > channel >= 0:
        adc_baseline = int(record.baseline[channel])
    adc_res = None
    if getattr(record, "adc_res", None) is not None and len(record.adc_res) > channel >= 0:
        adc_res = int(record.adc_res[channel])

    channel_signal, qa_metadata = validate_ecg_metadata_and_signal(
        signal_matrix=record.p_signal,
        channel=channel,
        fs=float(record.fs),
        units_list=list(record.units),
        sig_names=list(record.sig_name),
        r_peaks_ref=r_peaks_ref,
        expected_fs=expected_fs,
        expected_unit="mV",
        max_abs_mv=20.0,
        require_annotations=True,
        adc_gain=adc_gain,
        adc_baseline=adc_baseline,
        adc_res=adc_res,
    )

    # Record beat symbol breakdown and header comments for traceability
    symbol_counts: dict[str, int] = {}
    for s in beat_syms:
        symbol_counts[s] = symbol_counts.get(s, 0) + 1

    qa_metadata.update(
        {
            "record_id": rec_str,
            "database": pn_dir,
            "total_record_samples": total_samples,
            "all_annotations_in_window": int(len(beat_samples) + non_beat_count),
            "non_beat_annotations_excluded": int(non_beat_count),
            "beat_symbol_distribution": symbol_counts,
            "header_comments": list(getattr(record, "comments", []) or []),
        }
    )

    return ECGRecordData(
        record_id=rec_str,
        channel=channel,
        channel_name=qa_metadata["channel_name"],
        fs=float(record.fs),
        units=qa_metadata["physical_unit"],
        signal=channel_signal,
        r_peaks_ref=r_peaks_ref,
        beat_symbols=beat_syms,
        all_annotation_count=int(len(beat_samples) + non_beat_count),
        non_beat_annotation_count=int(non_beat_count),
        metadata=qa_metadata,
    )
