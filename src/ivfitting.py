import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


MZI_FSR_1550_NM = 14.3
MZI_FSR_1310_NM = 9.87
WAVELENGTH_BAND_SPLIT_NM = 1450.0
REFERENCE_POLY_DEGREE = 3
MZI_FIT_FSR_BOUNDS_SCALE = (0.65, 1.35)


def read_array(data_string):
    if data_string is None or data_string.strip() == '':
        return np.array([])
    return np.array([float(x) for x in data_string.split(',') if x.strip()])


def local_name(tag):
    return tag.split('}')[-1]


def child_text(element, child_name):
    for child in element:
        if local_name(child.tag) == child_name:
            return child.text
    return None


def first_descendant(element, tag_name):
    for child in element.iter():
        if local_name(child.tag) == tag_name:
            return child
    return None


def all_descendants(element, tag_name):
    return [child for child in element.iter() if local_name(child.tag) == tag_name]


def output_dir(filename):
    source_path = Path(filename)
    try:
        relative_path = source_path.relative_to('data')
    except ValueError:
        parts = source_path.parts
        relative_path = Path(*parts[1:]) if parts and parts[0] == 'data' else source_path.name

    return Path('res') / 'png' / relative_path.with_suffix('')


def finish_plot(fig, filename, save, show, plot_name):
    if save == 'T':
        path = output_dir(filename) / f'{plot_name}.png'
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path)
    if show == 'T':
        plt.show()
    plt.close(fig)


def combine_result_images(filename, save, show):
    result_dir = output_dir(filename)
    image_names = [
        '01_transmission_spectra_as_measured.png',
        '02_reference_polynomial_fits.png',
        '03_two_step_flattened_spectra.png',
        '04_mzi_fit_minus_1v.png',
        '05_iv_characteristics.png',
        '06_iv_analysis.png',
    ]
    image_paths = [result_dir / name for name in image_names]
    if not all(path.exists() for path in image_paths):
        missing = [path.name for path in image_paths if not path.exists()]
        print(f'{filename}: combined result was skipped. Missing files: {missing}')
        return

    titles = [
        'Transmission Spectra',
        'Polynomial Fits',
        'Two-step Flattened',
        'MZI Fit (-1.0V)',
        'IV Characteristics',
        'IV Analysis',
    ]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    for ax, image_path, title in zip(axes.ravel(), image_paths, titles):
        ax.imshow(plt.imread(image_path))
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.axis('off')

    fig.suptitle(Path(filename).stem, fontsize=14, fontweight='bold')
    fig.tight_layout()
    if save == 'T':
        result_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(result_dir / '00_combined_results.png')
    if show == 'T':
        plt.show()
    plt.close(fig)


def r_squared(y_true, y_fit):
    ss_res = np.sum((y_true - y_fit) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot == 0:
        return 0.0
    return 1 - (ss_res / ss_tot)


def smooth(y, window=5):
    if len(y) < window or window < 2:
        return y
    kernel = np.ones(window) / window
    return np.convolve(y, kernel, mode='same')


def clean_xy(x_vals, y_vals):
    size = min(len(x_vals), len(y_vals))
    if size == 0:
        return np.array([]), np.array([])

    x_vals = np.asarray(x_vals[:size], dtype=float)
    y_vals = np.asarray(y_vals[:size], dtype=float)
    finite = np.isfinite(x_vals) & np.isfinite(y_vals)
    x_vals = x_vals[finite]
    y_vals = y_vals[finite]

    if len(x_vals) == 0:
        return x_vals, y_vals

    order = np.argsort(x_vals)
    return x_vals[order], y_vals[order]


def fit_centered_polynomial(x_vals, y_vals, degree=REFERENCE_POLY_DEGREE):
    x_vals, y_vals = clean_xy(x_vals, y_vals)
    if len(x_vals) < 2:
        return None

    fit_degree = min(degree, len(x_vals) - 1)
    center = float(np.mean(x_vals))
    poly = np.poly1d(np.polyfit(x_vals - center, y_vals, fit_degree))
    return lambda wavelengths: poly(np.asarray(wavelengths) - center)


def restrict_to_wavelength_overlap(ref_l, ref_il, l_vals, il_vals):
    ref_l, ref_il = clean_xy(ref_l, ref_il)
    l_vals, il_vals = clean_xy(l_vals, il_vals)
    if len(ref_l) == 0 or len(l_vals) == 0:
        return ref_l, ref_il, l_vals, il_vals

    low = max(float(np.min(ref_l)), float(np.min(l_vals)))
    high = min(float(np.max(ref_l)), float(np.max(l_vals)))
    if high <= low:
        return np.array([]), np.array([]), np.array([]), np.array([])

    ref_mask = (ref_l >= low) & (ref_l <= high)
    data_mask = (l_vals >= low) & (l_vals <= high)
    return ref_l[ref_mask], ref_il[ref_mask], l_vals[data_mask], il_vals[data_mask]


def adaptive_peak_prominence(y_vals):
    if len(y_vals) < 3:
        return 0.5

    spread = float(np.nanpercentile(y_vals, 95) - np.nanpercentile(y_vals, 5))
    return max(0.1, min(1.0, spread * 0.05))


def robust_linear_peak_fit(peak_x, peak_y):
    if len(peak_x) < 2:
        return None, peak_x, peak_y

    if len(peak_x) < 4:
        return np.poly1d(np.polyfit(peak_x, peak_y, 1)), peak_x, peak_y

    rough_poly = np.poly1d(np.polyfit(peak_x, peak_y, 1))
    residual = peak_y - rough_poly(peak_x)
    median_residual = float(np.median(residual))
    mad = float(np.median(np.abs(residual - median_residual)))

    if mad == 0:
        keep = np.ones(len(peak_x), dtype=bool)
    else:
        keep = np.abs(residual - median_residual) <= 3.0 * 1.4826 * mad

    if np.count_nonzero(keep) < 2:
        keep = np.ones(len(peak_x), dtype=bool)

    return np.poly1d(np.polyfit(peak_x[keep], peak_y[keep], 1)), peak_x[keep], peak_y[keep]


def detect_envelope_peaks(l_vals, y_vals, expected_fsr_nm=None, distance_nm=None, prominence=None):
    if len(l_vals) < 3:
        return np.array([]), np.array([])

    y_smooth = smooth(y_vals)
    avg_step = abs(float(np.mean(np.diff(l_vals))))
    if distance_nm is None:
        distance_nm = expected_fsr_nm * 0.65 if expected_fsr_nm is not None else 10.0
    if prominence is None:
        prominence = adaptive_peak_prominence(y_vals)
    distance_pts = max(1, int(distance_nm / avg_step)) if avg_step > 0 else 1

    peak_idx, _ = find_peaks(y_smooth, distance=distance_pts, prominence=prominence)
    if len(peak_idx) == 0:
        return np.array([]), np.array([])

    peak_x = l_vals[peak_idx]
    peak_y = y_vals[peak_idx]
    return peak_x, peak_y


def flatten_spectrum(l_vals, il_vals, poly_ref_model, expected_fsr_nm=None):
    l_vals, il_vals = clean_xy(l_vals, il_vals)
    if poly_ref_model is None or len(l_vals) == 0:
        return np.array([]), None, (np.array([]), np.array([]))

    first_flattened = il_vals - poly_ref_model(l_vals)
    peak_x, peak_y = detect_envelope_peaks(l_vals, first_flattened, expected_fsr_nm=expected_fsr_nm)

    if len(peak_x) >= 2:
        residual_poly, fit_peak_x, fit_peak_y = robust_linear_peak_fit(peak_x, peak_y)
        return first_flattened - residual_poly(l_vals), residual_poly, (fit_peak_x, fit_peak_y)

    return first_flattened, None, (peak_x, peak_y)


def flatten_with_reference(ref_l, ref_il, l_vals, il_vals):
    ref_fit_l, ref_fit_il, l_vals, il_vals = restrict_to_wavelength_overlap(ref_l, ref_il, l_vals, il_vals)
    if len(ref_fit_l) < 2 or len(l_vals) < 3:
        return None

    poly_ref_model = fit_centered_polynomial(ref_fit_l, ref_fit_il)
    expected_fsr_nm = fsr_for_wavelength_band(l_vals)
    flattened, residual_poly, peaks = flatten_spectrum(
        l_vals,
        il_vals,
        poly_ref_model,
        expected_fsr_nm=expected_fsr_nm,
    )
    return l_vals, flattened, residual_poly, peaks


def get_wavelength_data(sweep):
    l_data = read_array(child_text(sweep, 'L'))
    il_data = read_array(child_text(sweep, 'IL'))
    if len(l_data) == 0 or len(il_data) == 0:
        return None, None
    return l_data, il_data


def find_modulator(root, exact_name=None, align=False):
    for elem in root.iter():
        if local_name(elem.tag) != 'Modulator':
            continue
        name = elem.attrib.get('Name', '')
        if exact_name is not None and name == exact_name:
            return elem
        if align and name.endswith('_ALIGN'):
            return elem
    return None


def get_align_and_mzm_data(root):
    align = find_modulator(root, exact_name='DCM_LMZC_ALIGN')
    if align is None:
        align = find_modulator(root, align=True)

    mzm = find_modulator(root, exact_name='MZMCTE_LULAB_450_500')
    if mzm is None:
        for elem in root.iter():
            if local_name(elem.tag) != 'Modulator':
                continue
            sweeps = all_descendants(elem, 'WavelengthSweep')
            if any('DCBias' in sweep.attrib for sweep in sweeps):
                mzm = elem
                break

    ref_l, ref_il = np.array([]), np.array([])
    if align is not None:
        ref_sweep = first_descendant(align, 'WavelengthSweep')
        if ref_sweep is not None:
            ref_l, ref_il = get_wavelength_data(ref_sweep)

    mzm_data = []
    if mzm is not None:
        for sweep in all_descendants(mzm, 'WavelengthSweep'):
            l_data, il_data = get_wavelength_data(sweep)
            if l_data is not None:
                bias = sweep.attrib.get('DCBias', '0.0')
                mzm_data.append((bias, l_data, il_data))

    return ref_l, ref_il, mzm_data


def plot_transmission_spectra(root, filename, save, show):
    fig = plt.figure(figsize=(8, 5))
    ref_l = None
    ref_il = None

    for elem in root.iter():
        if local_name(elem.tag) != 'WavelengthSweep':
            continue

        l_data, il_data = get_wavelength_data(elem)
        if l_data is None:
            continue

        if 'DCBias' not in elem.attrib:
            ref_l = l_data
            ref_il = il_data
        else:
            plt.plot(l_data, il_data, label=f"{elem.attrib.get('DCBias')} V")

    if ref_l is not None and ref_il is not None:
        plt.plot(ref_l, ref_il, color='pink', label='Reference')
        poly_model = np.poly1d(np.polyfit(ref_l, ref_il, 6))
        plt.plot(ref_l, poly_model(ref_l), 'k--', label='ref polynomial fit')

    plt.title('Transmission spectra - as measured', fontsize=14, fontweight='bold')
    plt.xlabel('Wavelength [nm]', fontsize=12)
    plt.ylabel('Measured transmission [dB]', fontsize=12)
    plt.grid(False)
    plt.legend(loc='lower left', ncol=2, fontsize=9)
    plt.tight_layout()
    finish_plot(fig, filename, save, show, '01_transmission_spectra_as_measured')


def plot_polynomial_fits(root, filename, save, show):
    ref_l = None
    ref_il = None

    for elem in root.iter():
        if local_name(elem.tag) != 'WavelengthSweep':
            continue

        dc_bias = elem.attrib.get('DCBias', 'Reference')
        l_data, il_data = get_wavelength_data(elem)
        if l_data is not None and dc_bias == '0.0':
            ref_l = l_data
            ref_il = il_data

    if ref_l is None or ref_il is None:
        print(f'{filename}: 0.0V reference data was not found.')
        return

    mean_l = np.mean(ref_l)
    l_centered = ref_l - mean_l

    fig = plt.figure(figsize=(9, 7))
    plt.plot(ref_l, ref_il, label='Raw Reference (0.0 V)', color='black', alpha=0.6, linewidth=2)

    for degree, color in zip([2, 3, 4, 5, 6], ['blue', 'red', 'green', 'orange', 'purple']):
        poly_model = np.poly1d(np.polyfit(l_centered, ref_il, degree))
        fit_il = poly_model(l_centered)
        r2 = r_squared(ref_il, fit_il)
        print(f'{filename}: {degree} degree polynomial fit R^2 = {r2:.6f}')
        plt.plot(
            ref_l,
            fit_il,
            label=f'{degree}th Degree Fit (R^2={r2:.4f})',
            color=color,
            linestyle='--',
            linewidth=1.5,
        )

    plt.title('Raw Reference Data vs Polynomial Fits (2nd~6th)', fontsize=14, fontweight='bold')
    plt.xlabel('Wavelength [nm]', fontsize=12)
    plt.ylabel('Transmission [dB]', fontsize=12)
    plt.legend(loc='lower left', fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    finish_plot(fig, filename, save, show, '02_reference_polynomial_fits')


def two_step_flatten(root):
    ref_l, ref_il, mzm_data = get_align_and_mzm_data(root)
    if len(ref_l) == 0 or len(ref_il) == 0 or not mzm_data:
        return ref_l, ref_il, [], None

    ref_l, ref_il = clean_xy(ref_l, ref_il)
    poly_ref_model = fit_centered_polynomial(ref_l, ref_il)
    if poly_ref_model is None:
        return ref_l, ref_il, [], None

    flattened_data = []

    for bias, l_vals, il_vals in mzm_data:
        flattened_result = flatten_with_reference(ref_l, ref_il, l_vals, il_vals)
        if flattened_result is None:
            continue
        l_vals, flattened, residual_poly, peaks = flattened_result
        flattened_data.append((bias, l_vals, flattened, residual_poly, peaks))

    return ref_l, ref_il, flattened_data, poly_ref_model


def plot_two_step_flattened(root, filename, save, show):
    ref_l, ref_il, flattened_data, poly_ref_model = two_step_flatten(root)
    if poly_ref_model is None:
        print(f'{filename}: data for two-step flattening was not found.')
        return

    fig = plt.figure(figsize=(12, 8))
    flattened_values = []

    if len(ref_l) > 0:
        ref_final, _, (ref_peak_x, ref_peak_y) = flatten_spectrum(
            ref_l,
            ref_il,
            poly_ref_model,
            expected_fsr_nm=fsr_for_wavelength_band(ref_l),
        )
        flattened_values.append(ref_final)
        plt.plot(ref_l, ref_final, color='black', linewidth=3, label='REF (2-step Flattened)', zorder=10)
        plt.scatter(ref_peak_x, ref_peak_y, color='black', marker='x', s=40, zorder=11)

    for bias, l_vals, two_step_flattened, residual_poly, (peak_x, peak_y) in flattened_data:
        flattened_values.append(two_step_flattened)
        line, = plt.plot(l_vals, two_step_flattened, label=f'Bias {bias}V', linewidth=1, alpha=0.8)
        color = line.get_color()
        if residual_poly is not None:
            plt.plot(l_vals, residual_poly(l_vals), '--', color=color, alpha=0.2, linewidth=0.8)
            plt.scatter(peak_x, peak_y, color=color, marker='o', s=18, alpha=0.6, zorder=5)

    plt.axhline(0, color='red', linewidth=0.5, linestyle='--')
    plt.title('Two-step Flattened Spectra (Band-aware Peak Detection)', fontsize=14, fontweight='bold')
    plt.xlabel('Wavelength [nm]', fontsize=12)
    plt.ylabel('Normalized Transmission [dB]', fontsize=12)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize='x-small', ncol=2)
    plt.grid(True, linestyle=':', alpha=0.6)
    if flattened_values:
        finite_value_sets = [values[np.isfinite(values)] for values in flattened_values if len(values) > 0]
        finite_value_sets = [values for values in finite_value_sets if len(values) > 0]
        if finite_value_sets:
            finite_values = np.concatenate(finite_value_sets)
            low = max(-80.0, float(np.percentile(finite_values, 1)) - 5.0)
            high = min(20.0, float(np.percentile(finite_values, 99)) + 5.0)
            plt.ylim(low, high)
    plt.tight_layout()
    finish_plot(fig, filename, save, show, '03_two_step_flattened_spectra')


def mzi_model(lam, a_value, b_value, lam_0, fsr, phi):
    return a_value + b_value * (np.cos(np.pi * (lam - lam_0) / fsr + phi)) ** 2


def fit_mzi(l_vals, flattened_db):
    l_vals, flattened_db = clean_xy(l_vals, flattened_db)
    if len(l_vals) < 5:
        print('MZI fitting failed: too few wavelength points.')
        return None, None, 0.0

    linear_spectrum = 10 ** (flattened_db / 10)
    max_val = np.max(linear_spectrum)
    if max_val <= 0 or not np.isfinite(max_val):
        print('MZI fitting failed: normalized spectrum maximum is 0.')
        return None, None, 0.0

    normalized_spectrum = linear_spectrum / max_val

    a_guess = np.min(normalized_spectrum)
    b_guess = np.max(normalized_spectrum) - a_guess
    wl0_guess = l_vals[np.argmax(normalized_spectrum)]
    fsr_guess = fsr_for_wavelength_band(l_vals) or 14.0
    phi_guess = 0.0
    fsr_lower = max(0.1, fsr_guess * MZI_FIT_FSR_BOUNDS_SCALE[0])
    fsr_upper = fsr_guess * MZI_FIT_FSR_BOUNDS_SCALE[1]

    initial_guess = [a_guess, b_guess, wl0_guess, fsr_guess, phi_guess]
    lower_bounds = [0.0, 0.0, float(np.min(l_vals)), fsr_lower, -np.pi]
    upper_bounds = [1.0, 2.0, float(np.max(l_vals)), fsr_upper, np.pi]

    try:
        popt, _ = curve_fit(
            mzi_model,
            l_vals,
            normalized_spectrum,
            p0=initial_guess,
            bounds=(lower_bounds, upper_bounds),
            maxfev=10000,
        )
    except Exception as exc:
        print(f'MZI fitting failed: {exc}')
        return None, normalized_spectrum, 0.0

    fit_curve = mzi_model(l_vals, *popt)
    return popt, normalized_spectrum, r_squared(normalized_spectrum, fit_curve)


def fsr_for_wavelength_band(wavelengths):
    if len(wavelengths) == 0:
        return None

    center_wavelength = float(np.mean(wavelengths))
    if center_wavelength >= WAVELENGTH_BAND_SPLIT_NM:
        return MZI_FSR_1550_NM
    return MZI_FSR_1310_NM


def parse_bias_value(bias):
    try:
        return float(bias)
    except (TypeError, ValueError):
        return None


def select_bias_data(flattened_data, target_bias=-1.0):
    numeric_candidates = []
    for item in flattened_data:
        bias_value = parse_bias_value(item[0])
        if bias_value is None:
            continue
        numeric_candidates.append((abs(bias_value - target_bias), bias_value, item))

    if not numeric_candidates:
        return None, None

    _, selected_bias, selected_item = min(numeric_candidates, key=lambda candidate: candidate[0])
    return selected_bias, selected_item


def plot_mzi_fit(root, filename, save, show):
    _, _, flattened_data, _ = two_step_flatten(root)
    selected_bias, selected_data = select_bias_data(flattened_data, target_bias=-1.0)

    if selected_data is None:
        print(f'{filename}: numeric bias flattened data was not found.')
        return

    bias, target_lam, target_flattened, *_ = selected_data
    popt, normalized_spectrum, r2 = fit_mzi(target_lam, target_flattened)
    if popt is None:
        print(f'{filename}: MZI fitting failed. Check initial values or bounds.')
        return

    fit_curve = mzi_model(target_lam, *popt)
    print(
        f'{filename}: MZI {selected_bias:.1f}V fit A={popt[0]:.4f}, B={popt[1]:.4f}, '
        f'lambda0={popt[2]:.4f} nm, FSR={popt[3]:.4f} nm, '
        f'phi={popt[4]:.4f}, R^2={r2:.4f}'
    )

    fig = plt.figure(figsize=(10, 6))
    plt.plot(
        target_lam,
        normalized_spectrum,
        color='C0',
        linestyle='-',
        alpha=0.75,
        label=f'Normalized Data ({selected_bias:.1f}V)',
    )
    plt.plot(target_lam, fit_curve, color='black', linestyle='--', label=f'MZI Fit (R^2 = {r2:.4f})', linewidth=2)
    plt.title(f'MZI Fitting - {selected_bias:.1f}V (Linear Scale)', fontsize=14, fontweight='bold')
    plt.xlabel('Wavelength [nm]', fontsize=12)
    plt.ylabel('Normalized Transmission [a.u.]', fontsize=12)
    plt.xlim(target_lam[0], target_lam[-1])
    plt.ylim(-0.05, 1.15)
    plt.legend(loc='upper right')
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()
    finish_plot(fig, filename, save, show, '04_mzi_fit_minus_1v')


def get_iv_data(root):
    for elem in root.iter():
        if local_name(elem.tag) == 'IVMeasurement':
            voltage = read_array(child_text(elem, 'Voltage'))
            current = read_array(child_text(elem, 'Current'))
            if len(voltage) > 0 and len(current) > 0:
                return voltage, current
    return None, None


def plot_iv_characteristics(root, filename, save, show):
    voltage, current = get_iv_data(root)
    if voltage is None:
        print(f'{filename}: IV measurement data was not found.')
        return

    fig = plt.figure(figsize=(8, 6))
    plt.semilogy(voltage, np.abs(current), 'bo', markersize=4, linewidth=1.5, label='|I|-V Curve')
    plt.title('IV Characteristics of Main Modulator', fontsize=14, fontweight='bold')
    plt.xlabel('Voltage [V]', fontsize=12)
    plt.ylabel('Absolute Current |I| [A]', fontsize=12)
    plt.grid(True, which='both', ls='--', alpha=0.5)
    plt.legend(loc='best')
    plt.tight_layout()
    finish_plot(fig, filename, save, show, '05_iv_characteristics')


def diode_model(voltage, i_s, n_value):
    thermal_voltage = 0.02585
    return i_s * (np.exp(voltage / (n_value * thermal_voltage)) - 1)


def plot_iv_analysis(root, filename, save, show):
    voltage, current = get_iv_data(root)
    if voltage is None:
        print(f'{filename}: IV measurement data was not found.')
        return

    abs_current = np.abs(current)
    rev_mask = voltage <= 0.25
    fwd_mask = voltage >= 0.5
    v_rev, i_rev = voltage[rev_mask], abs_current[rev_mask]
    v_fwd, i_fwd = voltage[fwd_mask], abs_current[fwd_mask]

    if len(v_rev) < 5 or len(v_fwd) < 2:
        print(f'{filename}: IV fitting range has too few data points.')
        return

    poly_model = np.poly1d(np.polyfit(v_rev, i_rev, 4))
    v_rev_fine = np.linspace(min(v_rev), max(v_rev), 100)
    fit_i_rev = poly_model(v_rev_fine)
    r2_rev = r_squared(i_rev, poly_model(v_rev))

    try:
        popt, _ = curve_fit(
            diode_model,
            v_fwd,
            i_fwd,
            p0=[1e-12, 1.5],
            bounds=([1e-16, 1.0], [1e-9, 5.0]),
            maxfev=100000,
        )
    except Exception as exc:
        print(f'{filename}: forward IV fitting failed: {exc}')
        return

    i_s_fit, n_fit = popt
    v_fwd_fine = np.linspace(min(v_fwd), max(v_fwd), 100)
    fit_i_fwd = diode_model(v_fwd_fine, i_s_fit, n_fit)
    r2_fwd = r_squared(i_fwd, diode_model(v_fwd, *popt))

    fig = plt.figure(figsize=(9, 7))
    plt.semilogy(voltage, abs_current, 'o', label='Measured IV', markersize=7)
    plt.semilogy(v_rev_fine, fit_i_rev, '-', color='tab:blue', linewidth=2.5, label='Reverse polynomial fit')
    plt.semilogy(v_fwd_fine, fit_i_fwd, '-', color='tab:red', linewidth=2.5, label='Forward diode fit')

    textstr = '\n'.join((
        f'Is = {i_s_fit:.3e} A',
        f'n = {n_fit:.3f}',
        f'R^2_fwd = {r2_fwd:.4f}',
        f'R^2_rev = {r2_rev:.4f}',
    ))
    plt.gca().text(
        0.05,
        0.95,
        textstr,
        transform=plt.gca().transAxes,
        fontsize=12,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
    )

    plt.title('IV analysis', fontsize=14)
    plt.xlabel('Voltage [V]', fontsize=12)
    plt.ylabel('Current [A]', fontsize=12)
    plt.grid(True, which='both', ls='-', alpha=0.3)
    plt.legend(loc='lower left')
    plt.tight_layout()
    finish_plot(fig, filename, save, show, '06_iv_analysis')


def IV(filename, save, show):
    tree = ET.parse(filename)
    root = tree.getroot()

    plot_transmission_spectra(root, filename, save, 'F')
    plot_polynomial_fits(root, filename, save, 'F')
    plot_two_step_flattened(root, filename, save, 'F')
    plot_mzi_fit(root, filename, save, 'F')
    plot_iv_characteristics(root, filename, save, 'F')
    plot_iv_analysis(root, filename, save, 'F')
    combine_result_images(filename, save, show)
