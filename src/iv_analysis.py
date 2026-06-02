from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit

from .config import THERMAL_VOLTAGE
from .xml_parser import r2_score


def diode_eq(voltage: np.ndarray, Is: float, n: float) -> np.ndarray:
    exponent = np.clip(voltage / (n * THERMAL_VOLTAGE), -700, 700)
    return Is * (np.exp(exponent) - 1.0)


def plot_iv_log(ax, iv: dict[str, np.ndarray]) -> None:
    if iv["V"].size:
        current_abs = np.abs(iv["I"])
        positive = current_abs > 0
        ax.semilogy(iv["V"][positive], current_abs[positive], "bo", ms=4)
    ax.set_title("IV-curve (Log Scale)")
    ax.set_xlabel("Voltage [V]")
    ax.set_ylabel("Current [A]")
    ax.grid(True, which="both", ls="--", alpha=0.5)


def plot_iv_analysis(ax, iv: dict[str, np.ndarray]) -> None:
    if not iv["V"].size:
        ax.set_title("IV analysis")
        ax.set_xlabel("Voltage [V]")
        ax.set_ylabel("Current [A]")
        ax.grid(True, which="both", ls="--", alpha=0.5)
        ax.legend(fontsize="small", loc="lower left")
        return

    voltage = iv["V"]
    current_abs = np.abs(iv["I"])
    positive = current_abs > 0
    ax.semilogy(voltage[positive], current_abs[positive], "o",
                color="tab:blue", ms=5, label="Measured IV")

    reverse = (voltage < 0.5) & positive
    r2_rev = float("nan")
    if np.count_nonzero(reverse) >= 4:
        v_rev_pts = voltage[reverse]
        log_i_rev = np.log10(current_abs[reverse])
        keep = np.ones(v_rev_pts.size, dtype=bool)
        for i in range(1, v_rev_pts.size - 1):
            if log_i_rev[i] < 0.5 * (log_i_rev[i - 1] + log_i_rev[i + 1]) - 1.0:
                keep[i] = False
        v_keep, log_keep = v_rev_pts[keep], log_i_rev[keep]
        deg = min(4, v_keep.size - 1)
        rev_poly = np.poly1d(np.polyfit(v_keep, log_keep, deg))
        v_rev = np.linspace(float(v_keep.min()), float(v_keep.max()), 200)
        ax.semilogy(v_rev, 10 ** rev_poly(v_rev), "-", color="tab:orange",
                    label="Reverse polynomial fit")
        r2_rev = r2_score(log_keep, rev_poly(v_keep))

    forward = (voltage >= 0.5) & positive
    diode_is = diode_n = r2_fwd = float("nan")
    turn_on = (np.count_nonzero(forward) >= 2 and
               (np.log10(current_abs[forward].max()) -
                np.log10(current_abs[forward].min())) > 1.0)
    if turn_on:
        try:
            popt, _ = curve_fit(
                diode_eq, voltage[forward], current_abs[forward],
                p0=[1e-15, 1.5], bounds=([1e-30, 0.5], [1e-3, 10.0]), maxfev=10000,
            )
            diode_is, diode_n = float(popt[0]), float(popt[1])
            v_fwd = np.linspace(float(voltage[forward].min()), float(voltage.max()), 100)
            ax.semilogy(v_fwd, diode_eq(v_fwd, *popt), "-", color="tab:green",
                        label="Forward diode fit")
            r2_fwd = r2_score(current_abs[forward], diode_eq(voltage[forward], *popt))
        except Exception as exc:
            ax.text(0.5, 0.2, f"Diode fit failed\n{exc}", transform=ax.transAxes,
                    ha="center", color="red", fontsize="small")

    def _fmt(value, spec):
        return "n/a" if not np.isfinite(value) else format(value, spec)

    note = "" if turn_on else "\n(flat IV: no turn-on)"
    stats = (f"Is = {_fmt(diode_is, '.3e')} A\nn = {_fmt(diode_n, '.3f')}\n"
             f"$R^2_{{fwd}}$ = {_fmt(r2_fwd, '.4f')}\n$R^2_{{rev}}$ = {_fmt(r2_rev, '.4f')}{note}")
    ax.text(0.03, 0.97, stats, transform=ax.transAxes, va="top", ha="left",
            fontsize=9, family="monospace",
            bbox=dict(boxstyle="square,pad=0.4", fc="white", ec="0.4", lw=0.8))

    ax.set_title("IV analysis")
    ax.set_xlabel("Voltage [V]")
    ax.set_ylabel("Current [A]")
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.legend(fontsize="small", loc="lower left")
