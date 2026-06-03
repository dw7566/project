from __future__ import annotations

import xml.etree.ElementTree as ET

import numpy as np

from .xml_parser import attr_any


O_BAND_RANGE_NM = (1260.0, 1360.0)
C_BAND_RANGE_NM = (1530.0, 1580.0)


def insertion_loss_band(*names: object) -> tuple[str, tuple[float, float] | None]:
    text = " ".join(str(name).upper() for name in names if name)
    if "LMZO" in text:
        return "O-band", O_BAND_RANGE_NM
    if "LMZC" in text:
        return "C-band", C_BAND_RANGE_NM
    return "Full band", None


def insertion_loss_band_from_root(root: ET.Element) -> tuple[str, tuple[float, float] | None]:
    test_site_info = root.find(".//TestSiteInfo")
    names: list[object] = [
        attr_any(test_site_info, "TestSite"),
        attr_any(test_site_info, "Maskset"),
    ]
    for device_info in root.findall(".//DeviceInfo"):
        names.append(device_info.get("Name"))
    return insertion_loss_band(*names)


def positive_insertion_loss(il_db: np.ndarray) -> np.ndarray:
    finite = il_db[np.isfinite(il_db)]
    if finite.size and float(np.nanmean(finite)) <= 0.0:
        return -il_db
    return il_db


def band_mask(
    wavelength: np.ndarray,
    band_range: tuple[float, float] | None,
) -> np.ndarray:
    mask = np.isfinite(wavelength)
    if band_range is None:
        return mask
    low, high = band_range
    selected = mask & (wavelength >= low) & (wavelength <= high)
    return selected if np.count_nonzero(selected) else mask


def insertion_loss_db(
    wavelength: list[float],
    il: list[float],
    *names: object,
) -> float | None:
    count = min(len(wavelength), len(il))
    if count == 0:
        return None

    wavelength_arr = np.asarray(wavelength[:count], dtype=float)
    loss_arr = positive_insertion_loss(np.asarray(il[:count], dtype=float))
    _, band_range = insertion_loss_band(*names)
    keep = band_mask(wavelength_arr, band_range) & np.isfinite(loss_arr)
    if not np.any(keep):
        return None
    return float(np.nanmean(loss_arr[keep]))


def plot_insertion_loss_panel(
    ax,
    root: ET.Element,
    sweeps: list[dict[str, object]],
) -> None:
    band_name, band_range = insertion_loss_band_from_root(root)
    plotted = False

    for index, sweep in enumerate(sweeps):
        wavelength = sweep["L"]
        il = sweep["IL"]
        if not isinstance(wavelength, np.ndarray) or not isinstance(il, np.ndarray):
            continue

        count = min(wavelength.size, il.size)
        if count < 2:
            continue

        wavelength = wavelength[:count]
        loss = positive_insertion_loss(il[:count])
        keep = band_mask(wavelength, band_range) & np.isfinite(loss)
        if np.count_nonzero(keep) < 2:
            continue

        bias = str(sweep.get("Bias", ""))
        label = f"Reference ({bias}V)" if index == len(sweeps) - 1 else f"Bias {bias}V"
        ax.plot(wavelength[keep], loss[keep], linewidth=1.0, label=label)
        plotted = True

    if not plotted:
        ax.set_axis_off()
        ax.text(
            0.5,
            0.5,
            "Insertion loss\nNo wavelength sweep data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="red",
        )
        return

    if band_range is not None:
        ax.set_xlim(*band_range)
    ax.set_title(f"Insertion loss vs wavelength ({band_name})")
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Insertion loss [dB]")
    ax.grid(True, ls="--", alpha=0.35)
    ax.legend(ncol=3, fontsize="x-small", loc="best")
