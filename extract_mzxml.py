from pyopenms import MSExperiment, MzXMLFile
import os
import numpy as np
import pandas as pd
import traceback
from matplotlib.backends.backend_pdf import PdfPages
from extract_mzxml import *  # Ensure this contains necessary extraction functions
from plot_scannumber import plot_spectrum
from pyteomics import mzxml
from pyopenms import MSExperiment, MzXMLFile



def get_baseline(inten):
    import numpy as np
    if not isinstance(inten, np.ndarray):
        inten = np.array(inten)
    # Avoid log(0) by adding a small constant
    log_inten = np.log10(inten + 1e-9)
    t = np.arange(np.min(log_inten), np.max(log_inten) + 0.08, 0.08)
    n, xout = np.histogram(log_inten, bins=t)
    idx = np.argmax(n)
    baseline = 10 ** xout[idx]
    return baseline


def extract_ms1_scans(mzxml_file):

    scan_data = []

    exp = MSExperiment()
    MzXMLFile().load(mzxml_file, exp)
    ms1_spectra = [s for s in exp.getSpectra() if s.getMSLevel() == 1]

    for spectrum in ms1_spectra:
        scan_number = spectrum.getNativeID()
        mz, intensity = spectrum.get_peaks()
        mz = np.array(mz, dtype=float)
        intensity = np.array(intensity, dtype=float)

        if len(intensity) == 0:
            continue

        # Apply baseline correction
        baseline = get_baseline(intensity)
        corrected_intensity = intensity - baseline
        corrected_intensity[corrected_intensity < 0] = 0

        # Calculate relative abundance
        max_intensity = np.max(corrected_intensity)
        relative_abundance = (corrected_intensity / max_intensity) * 100 if max_intensity != 0 else np.zeros_like(
            corrected_intensity)

        scan_data.append((scan_number, (mz, relative_abundance)))

    return scan_data


def save_scan_data(scan_data, mzxml_file):

    output_dir = os.path.dirname(mzxml_file)
    output_file = os.path.join(output_dir, "ms1_scan_data.csv")

    # Convert scan_data into a structured format for CSV
    formatted_data = []
    for scan_number, (mz_values, relative_abundances) in scan_data:
        # Ensure scan_number is stored as an integer
        try:
            scan_number = int(scan_number)
        except ValueError:
            scan_number = int(''.join(filter(str.isdigit, str(scan_number))))

        formatted_data.append({
            "Scan_Number": scan_number,
            "m/z Values": np.array(mz_values, dtype=np.float64).tolist(),  # Store as float array
            "Relative Abundances": np.array(relative_abundances, dtype=np.float64).tolist()  # Store as float array
        })

    df = pd.DataFrame(formatted_data)
    df.to_csv(output_file, index=False, float_format='%.6f')  # Ensure numeric precision
    print(f"MS1 scan data saved to: {output_file}")
