import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import correlate
from scipy.ndimage import gaussian_filter1d
import numpy as np
from scipy.signal import correlate
from scipy.ndimage import gaussian_filter1d

def calculate_similarity_score(theoretical_mz, theoretical_abundance, experimental_mz, experimental_abundance):

    # Interpolate to align the m/z values
    common_mz = np.union1d(theoretical_mz, experimental_mz)
    theoretical_interp = np.interp(common_mz, theoretical_mz, theoretical_abundance, left=0, right=0)
    experimental_interp = np.interp(common_mz, experimental_mz, experimental_abundance, left=0, right=0)

    # Normalize the distributions
    theoretical_interp /= np.linalg.norm(theoretical_interp)
    experimental_interp /= np.linalg.norm(experimental_interp)

    # Calculate cosine similarity
    cosine_similarity = np.dot(theoretical_interp, experimental_interp)
    return cosine_similarity

def calculate_peak_pattern_similarity(theoretical_mz, theoretical_abundance, experimental_mz, experimental_abundance):

    # Smooth the distributions to reduce noise
    theoretical_abundance_smoothed = gaussian_filter1d(theoretical_abundance, sigma=1)
    experimental_abundance_smoothed = gaussian_filter1d(experimental_abundance, sigma=1)

    # Normalize the distributions
    theoretical_abundance_smoothed /= np.max(theoretical_abundance_smoothed)
    experimental_abundance_smoothed /= np.max(experimental_abundance_smoothed)

    # Calculate cross-correlation
    correlation = correlate(theoretical_abundance_smoothed, experimental_abundance_smoothed, mode='same')
    peak_pattern_similarity = np.max(correlation) / np.sqrt(np.sum(theoretical_abundance_smoothed**2) * np.sum(experimental_abundance_smoothed**2))

    return peak_pattern_similarity

def calculate_rms_difference(theoretical_mz, theoretical_abundance, experimental_mz, experimental_abundance):

    # Interpolate to align the m/z values
    common_mz = np.union1d(theoretical_mz, experimental_mz)
    theoretical_interp = np.interp(common_mz, theoretical_mz, theoretical_abundance, left=0, right=0)
    experimental_interp = np.interp(common_mz, experimental_mz, experimental_abundance, left=0, right=0)

    # Calculate RMS difference
    rms_difference = np.sqrt(np.mean((theoretical_interp - experimental_interp)**2))
    return rms_difference

def calculate_peak_pattern_similarity(theoretical_mz, theoretical_abundance, experimental_mz, experimental_abundance):

    # Smooth the distributions to reduce noise
    theoretical_abundance_smoothed = gaussian_filter1d(theoretical_abundance, sigma=1)
    experimental_abundance_smoothed = gaussian_filter1d(experimental_abundance, sigma=1)

    # Normalize the distributions
    theoretical_abundance_smoothed /= np.max(theoretical_abundance_smoothed)
    experimental_abundance_smoothed /= np.max(experimental_abundance_smoothed)

    # Calculate cross-correlation
    correlation = correlate(theoretical_abundance_smoothed, experimental_abundance_smoothed, mode='same')
    similarity = np.max(correlation) / np.sqrt(np.sum(theoretical_abundance_smoothed**2) * np.sum(experimental_abundance_smoothed**2))

    return similarity

def calculate_peak_range(mz_values, abundances, threshold=0.05):

    if len(mz_values) == 0:
        return None, None

    # Find the index of the maximum abundance
    max_abundance_index = np.argmax(abundances)
    max_abundance = abundances[max_abundance_index]

    # Find the start of the peak (where abundance rises above the threshold)
    peak_start = mz_values[0]
    for i in range(max_abundance_index, -1, -1):
        if abundances[i] < max_abundance * threshold:
            peak_start = mz_values[i]
            break

    # Find the end of the peak (where abundance drops below the threshold)
    peak_end = mz_values[-1]
    for i in range(max_abundance_index, len(abundances)):
        if abundances[i] < max_abundance * threshold:
            peak_end = mz_values[i]
            break

    return peak_start, peak_end

def find_densest_region(mz_values, abundances, window_size=20):

    max_density = 0
    densest_start = mz_values[0]
    densest_end = mz_values[-1]

    for i in range(len(mz_values) - window_size):
        current_density = np.sum(abundances[i:i+window_size])
        if current_density > max_density:
            max_density = current_density
            densest_start = mz_values[i]
            densest_end = mz_values[i + window_size]

    return densest_start, densest_end


def save_stem_plot_to_pdf(mz_values, abundances, peptide_sequence, scan_number, pdf, zoom=False):

    plt.figure(figsize=(10, 6))

    # Stem plot with red peaks, NO circular markers
    plt.stem(mz_values, abundances, linefmt='r-', markerfmt=' ', basefmt=' ')

    # Labels and Title
    plt.xlabel("m/z")
    plt.ylabel("Relative Abundance")
    plt.title(f"Peptide: {peptide_sequence}")

    # Grid for better readability
    plt.grid(True)

    # Zoom in on the peak region if requested
    if zoom:
        densest_start, densest_end = find_densest_region(mz_values, abundances)
        if densest_start is not None and densest_end is not None:
            # Calculate the center of the densest region
            center = (densest_start + densest_end) / 2

            # Expand the window size further (e.g., 2.0 for 2x expansion)
            window_size = (densest_end - densest_start) * 5  # Increase multiplier for larger window
            plt.xlim(center - window_size, center + window_size)
            plt.title(f"Peptide: {peptide_sequence} (Zoomed to Densest Region)")

    # Force x-axis to show only integer values
    x_min, x_max = plt.xlim()  # Get current x-axis limits
    plt.xticks(np.arange(int(x_min), int(x_max) + 1, step=1))  # Set ticks to integer steps

    # Save figure to the interactive PDF
    pdf.savefig()
    plt.close()

    print(f"Plot for Scan {scan_number} added to PDF.")