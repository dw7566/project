from __future__ import annotations

import xml.etree.ElementTree as ET

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import find_peaks

from .xml_parser import attr_any


def mzi_model(wavelength: np.ndarray, A: float, B: float, wl0: float,
              FSR: float, phi: float, slope: float) -> np.ndarray:
    x = wavelength - wavelength.mean()
    return A + slope * x + B * np.cos(np.pi * (wavelength - wl0) / FSR + phi) ** 2


def measure_fsr(wavelength: np.ndarray, transmission_db: np.ndarray, fallback: float) -> float:
    step = float(np.median(np.diff(wavelength))) if wavelength.size > 2 else 0.01
    min_distance = max(1, int(round(2.0 / max(abs(step), 1e-6))))
    for prominence in (1.0, 0.5, 0.25):
        notches, _ = find_peaks(-transmission_db, prominence=prominence, distance=min_distance)
        if len(notches) >= 2:
            spacing = np.diff(wavelength[notches])
            spacing = spacing[np.isfinite(spacing) & (spacing > 0)]
            if spacing.size:
                return float(np.median(spacing))
    return fallback


def crest_points(wavelength: np.ndarray, values_db: np.ndarray,
                 fsr: float) -> tuple[np.ndarray, np.ndarray]:
    step = float(np.median(np.diff(wavelength))) if wavelength.size > 2 else 0.01
    distance = max(1, int(round(0.6 * fsr / max(abs(step), 1e-6))))
    span = float(np.percentile(values_db, 95) - np.percentile(values_db, 5))
    prominence = float(np.clip(0.15 * span, 2.0, 8.0))
    peaks, _ = find_peaks(values_db, distance=distance, prominence=prominence)
    xs = list(wavelength[peaks])
    ys = list(values_db[peaks])
    for lo, hi in ((wavelength.min(), wavelength.min() + fsr),
                   (wavelength.max() - fsr, wavelength.max())):
        mask = (wavelength >= lo) & (wavelength <= hi)
        if np.any(mask):
            idx = int(np.argmax(values_db[mask]))
            xs.append(float(wavelength[mask][idx]))
            ys.append(float(values_db[mask][idx]))
    xs_arr = np.asarray(xs)
    ys_arr = np.asarray(ys)
    order = np.argsort(xs_arr)
    xs_arr, ys_arr = xs_arr[order], ys_arr[order]
    _, unique_idx = np.unique(np.round(xs_arr, 3), return_index=True)
    xs_arr, ys_arr = xs_arr[unique_idx], ys_arr[unique_idx]
    for _ in range(4):
        if xs_arr.size < 4:
            break
        trend = np.poly1d(np.polyfit(xs_arr, ys_arr, 1))(xs_arr)
        keep = (ys_arr - trend) > -2.5
        if keep.all():
            break
        xs_arr, ys_arr = xs_arr[keep], ys_arr[keep]
    return xs_arr, ys_arr


def top_envelope(wavelength: np.ndarray, values_db: np.ndarray,
                 fsr: float, degree: int = 2) -> np.ndarray:
    xs, ys = crest_points(wavelength, values_db, fsr)
    if xs.size >= 2:
        envelope = PchipInterpolator(xs, ys, extrapolate=True)(wavelength)
        envelope[wavelength < xs[0]] = ys[0]
        envelope[wavelength > xs[-1]] = ys[-1]
        return envelope
    n_windows = int(np.clip(round((wavelength.max() - wavelength.min()) / fsr), 3, 12))
    edges = np.linspace(float(wavelength.min()), float(wavelength.max()), n_windows + 1)
    xs_l: list[float] = []
    ys_l: list[float] = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (wavelength >= lo) & (wavelength <= hi)
        if np.any(mask):
            idx = int(np.argmax(values_db[mask]))
            xs_l.append(float(wavelength[mask][idx]))
            ys_l.append(float(values_db[mask][idx]))
    if len(xs_l) < 2:
        return np.zeros_like(wavelength)
    degree = min(degree, len(xs_l) - 1)
    return np.poly1d(np.polyfit(xs_l, ys_l, degree))(wavelength)


def flatten_to_envelope(wavelength: np.ndarray, values_db: np.ndarray,
                        fsr: float | None = None, degree: int = 2) -> np.ndarray:
    if wavelength.size < 4:
        return values_db
    band = float(wavelength.max() - wavelength.min())
    if not fsr or fsr <= 0:
        fsr = band / 6.0
    return values_db - top_envelope(wavelength, values_db, fsr, degree)


def device_fsr_fallback(root: ET.Element) -> float:
    site = attr_any(root.find(".//TestSiteInfo"), "TestSite").upper()
    if "LMZO" in site:
        return 10.0
    if "LMZC" in site:
        return 14.0
    return 12.0
