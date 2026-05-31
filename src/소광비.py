import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt

# ── 설정 ──────────────────────────────────────────────────
XML_PATH = r"C:\PythonProject\project\data\HY202103\D07\20190715_190855\HY202103_D07_(0,0)_LION1_DCM_LMZC.xml"
# 분석할 파장 범위 (nm). None 이면 전체 사용
LAMBDA_MIN = 1530.0
LAMBDA_MAX = 1580.0


# ──────────────────────────────────────────────────────────


def parse_wavelength_sweeps(xml_path: str) -> list[dict]:
    """XML에서 WavelengthSweep 데이터를 파싱합니다."""
    tree = ET.parse(xml_path)
    root = tree.getroot()

    sweeps = []
    for ws in root.iter("WavelengthSweep"):
        dc_bias = float(ws.get("DCBias", "0"))
        l_elem = ws.find("L")
        il_elem = ws.find("IL")

        if l_elem is None or il_elem is None:
            continue

        wavelengths = np.array([float(v) for v in l_elem.text.split(",")])
        insertion_loss = np.array([float(v) for v in il_elem.text.split(",")])
        sweeps.append({"dc_bias": dc_bias, "wavelengths": wavelengths, "IL_dBm": insertion_loss})

    return sweeps


def compute_extinction_ratio(wavelengths: np.ndarray,
                             il_dbm: np.ndarray,
                             lam_min: float | None = None,
                             lam_max: float | None = None) -> float:
    """
    Extinction Ratio (ER) = max(IL) - min(IL)  [단위: dB]

    파라미터
    --------
    wavelengths : 파장 배열 (nm)
    il_dbm      : 삽입 손실 배열 (dBm)
    lam_min/max : 분석 파장 범위 (None 이면 전체)

    반환
    ----
    er_dB : extinction ratio (dB, 양수)
    """
    mask = np.ones(len(wavelengths), dtype=bool)
    if lam_min is not None:
        mask &= wavelengths >= lam_min
    if lam_max is not None:
        mask &= wavelengths <= lam_max

    il_range = il_dbm[mask]
    if len(il_range) == 0:
        return float("nan")

    er = il_range.max() - il_range.min()  # dB (양수)
    return float(er)


def main():
    sweeps = parse_wavelength_sweeps(XML_PATH)

    print(f"{'DCBias (V)':>12}  {'ER (dB)':>10}  {'IL_max (dBm)':>14}  {'IL_min (dBm)':>14}")
    print("-" * 60)

    results = []
    for sw in sweeps:
        er = compute_extinction_ratio(
            sw["wavelengths"], sw["IL_dBm"], LAMBDA_MIN, LAMBDA_MAX
        )
        il_max = sw["IL_dBm"].max()
        il_min = sw["IL_dBm"].min()
        results.append((sw["dc_bias"], er, il_max, il_min))
        print(f"{sw['dc_bias']:>12.1f}  {er:>10.4f}  {il_max:>14.4f}  {il_min:>14.4f}")

    # ── 그래프 ──────────────────────────────────────────────
    biases = [r[0] for r in results]
    ers = [r[1] for r in results]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # (1) IL 스펙트럼 (각 bias)
    ax1 = axes[0]
    for sw in sweeps:
        mask = (sw["wavelengths"] >= (LAMBDA_MIN or -np.inf)) & \
               (sw["wavelengths"] <= (LAMBDA_MAX or np.inf))
        ax1.plot(sw["wavelengths"][mask], sw["IL_dBm"][mask],
                 label=f'{sw["dc_bias"]} V')
    ax1.set_xlabel("Wavelength (nm)")
    ax1.set_ylabel("Insertion Loss (dBm)")
    ax1.set_title("Transmission Spectrum by DC Bias")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # (2) ER vs Bias
    ax2 = axes[1]
    ax2.bar([str(b) for b in biases], ers, color="steelblue", edgecolor="navy")
    ax2.set_xlabel("DC Bias (V)")
    ax2.set_ylabel("Extinction Ratio (dB)")
    ax2.set_title("Extinction Ratio vs DC Bias")
    ax2.grid(axis="y", alpha=0.3)
    for i, er in enumerate(ers):
        ax2.text(i, er + 0.1, f"{er:.2f}", ha="center", va="bottom", fontsize=9)

    plt.tight_layout()
    plt.savefig("extinction_ratio.png", dpi=150)
    plt.show()
    print("\n그래프가 extinction_ratio.png 로 저장되었습니다.")


if __name__ == "__main__":