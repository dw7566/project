from __future__ import annotations

import xml.etree.ElementTree as ET

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import argrelextrema

from .xml_parser import attr_any


# ── 밴드 정의 ──
O_BAND_RANGE_NM = (1260.0, 1360.0)
C_BAND_RANGE_NM = (1530.0, 1580.0)
O_BAND_CENTER_NM = 1310.0
C_BAND_CENTER_NM = 1550.0


def insertion_loss_band(*names: object) -> tuple[str, tuple[float, float] | None, float | None]:
    """디바이스 이름으로 O/C-band 판별. (밴드명, 범위, 중심파장) 반환."""
    text = " ".join(str(name).upper() for name in names if name)
    if "LMZO" in text:
        return "O-band", O_BAND_RANGE_NM, O_BAND_CENTER_NM
    if "LMZC" in text:
        return "C-band", C_BAND_RANGE_NM, C_BAND_CENTER_NM
    return "Full band", None, None


def insertion_loss_band_from_root(root: ET.Element) -> tuple[str, tuple[float, float] | None, float | None]:
    test_site_info = root.find(".//TestSiteInfo")
    names: list[object] = [
        attr_any(test_site_info, "TestSite"),
        attr_any(test_site_info, "Maskset"),
    ]
    for device_info in root.findall(".//DeviceInfo"):
        names.append(device_info.get("Name"))
    return insertion_loss_band(*names)


def _band_mask(wavelength: np.ndarray, band_range: tuple[float, float] | None) -> np.ndarray:
    mask = np.isfinite(wavelength)
    if band_range is None:
        return mask
    low, high = band_range
    selected = mask & (wavelength >= low) & (wavelength <= high)
    return selected if np.count_nonzero(selected) else mask


def _peak_envelope(wavelength: np.ndarray, il_db: np.ndarray, order: int = 30) -> np.ndarray:
    """
    Local maxima를 PCHIP 보간하여 smooth peak envelope를 생성한다.
    MZM fringe의 peak(constructive interference) 레벨만 추적.
    """
    maxima_idx = argrelextrema(il_db, np.greater, order=order)[0]

    # 양 끝단 보정: 첫/마지막 FSR 구간에서 최대값 추가
    if maxima_idx.size >= 2:
        first_gap = maxima_idx[0]
        if first_gap > order:
            edge_idx = int(np.argmax(il_db[:first_gap]))
            maxima_idx = np.insert(maxima_idx, 0, edge_idx)
        last_gap = len(il_db) - 1 - maxima_idx[-1]
        if last_gap > order:
            edge_idx = maxima_idx[-1] + 1 + int(np.argmax(il_db[maxima_idx[-1] + 1:]))
            maxima_idx = np.append(maxima_idx, edge_idx)

    if maxima_idx.size < 2:
        return il_db.copy()

    # 중복 제거 + 정렬
    maxima_idx = np.unique(maxima_idx)
    peak_wl = wavelength[maxima_idx]
    peak_il = il_db[maxima_idx]

    # PCHIP 보간 (monotone, overshoot 없음)
    envelope = PchipInterpolator(peak_wl, peak_il, extrapolate=True)(wavelength)
    # 끝단 flat extrapolation
    envelope[wavelength < peak_wl[0]] = peak_il[0]
    envelope[wavelength > peak_wl[-1]] = peak_il[-1]

    return envelope


def _extract_il_from_envelope(
    wavelength: np.ndarray,
    il_db: np.ndarray,
    band_range: tuple[float, float] | None,
    center_nm: float | None,
) -> dict[str, float]:
    """
    Envelope 기반 IL 추출.
    Returns: {il_mean, il_at_center, il_min, il_max}  (모두 dBm 단위)
    """
    mask = _band_mask(wavelength, band_range)
    wl_band = wavelength[mask]
    il_band = il_db[mask]

    if wl_band.size < 4:
        return {
            "il_mean": float("nan"),
            "il_at_center": float("nan"),
            "il_min": float("nan"),
            "il_max": float("nan"),
        }

    envelope = _peak_envelope(wl_band, il_band)

    il_mean = float(np.nanmean(envelope))
    il_min = float(np.nanmin(envelope))
    il_max = float(np.nanmax(envelope))

    if center_nm is not None and wl_band[0] <= center_nm <= wl_band[-1]:
        il_at_center = float(np.interp(center_nm, wl_band, envelope))
    else:
        il_at_center = il_mean

    return {
        "il_mean": il_mean,
        "il_at_center": il_at_center,
        "il_min": il_min,
        "il_max": il_max,
    }


def insertion_loss_db(
    wavelength: list[float],
    il: list[float],
    *names: object,
) -> float | None:
    """
    Envelope 기반 삽입손실 (dBm).
    밴드 중심파장에서의 envelope 값을 반환한다.
    """
    count = min(len(wavelength), len(il))
    if count < 4:
        return None

    wl_arr = np.asarray(wavelength[:count], dtype=float)
    il_arr = np.asarray(il[:count], dtype=float)
    _, band_range, center_nm = insertion_loss_band(*names)

    result = _extract_il_from_envelope(wl_arr, il_arr, band_range, center_nm)
    value = result["il_at_center"]
    return value if np.isfinite(value) else None


def plot_insertion_loss_panel(
    ax,
    root: ET.Element,
    sweeps: list[dict[str, object]],
) -> None:
    """
    삽입손실 플롯: raw 스펙트럼 + peak envelope 오버레이.
    Reference sweep의 envelope을 강조 표시한다.
    """
    band_name, band_range, center_nm = insertion_loss_band_from_root(root)
    plotted = False

    for index, sweep in enumerate(sweeps):
        wavelength = sweep["L"]
        il = sweep["IL"]
        if not isinstance(wavelength, np.ndarray) or not isinstance(il, np.ndarray):
            continue

        count = min(wavelength.size, il.size)
        if count < 4:
            continue

        wavelength = wavelength[:count]
        il_db = il[:count]
        mask = _band_mask(wavelength, band_range)
        if np.count_nonzero(mask) < 4:
            continue

        wl_m = wavelength[mask]
        il_m = il_db[mask]

        bias = str(sweep.get("Bias", ""))
        is_ref = index == len(sweeps) - 1
        label = f"Reference ({bias}V)" if is_ref else f"Bias {bias}V"

        # raw 스펙트럼
        ax.plot(wl_m, il_m, linewidth=0.6, alpha=0.5, label=label)

        # Reference sweep에만 envelope 오버레이
        if is_ref:
            envelope = _peak_envelope(wl_m, il_m)
            ax.plot(wl_m, envelope, linewidth=2.0, color="black",
                    linestyle="--", label="Peak envelope", zorder=10)

            # 중심파장 IL 표시
            if center_nm is not None and wl_m[0] <= center_nm <= wl_m[-1]:
                il_center = float(np.interp(center_nm, wl_m, envelope))
                ax.axvline(center_nm, color="red", linewidth=0.8, alpha=0.5)
                ax.plot(center_nm, il_center, "ro", ms=8, zorder=11)
                ax.annotate(
                    f"IL = {il_center:.1f} dB\n@ {center_nm:.0f} nm",
                    xy=(center_nm, il_center),
                    xytext=(15, 15), textcoords="offset points",
                    fontsize=8, color="red",
                    arrowprops=dict(arrowstyle="->", color="red", lw=1.0),
                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="red",
                              boxstyle="round,pad=0.3"),
                )

        plotted = True

    if not plotted:
        ax.set_axis_off()
        ax.text(0.5, 0.5, "Insertion loss\nNo wavelength sweep data",
                transform=ax.transAxes, ha="center", va="center", color="red")
        return

    if band_range is not None:
        ax.set_xlim(*band_range)
    ax.set_title(f"Insertion loss ({band_name}) — envelope method")
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Transmission [dB]")
    ax.grid(True, ls="--", alpha=0.35)
    ax.legend(ncol=2, fontsize="x-small", loc="best")