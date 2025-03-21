import matplotlib.pyplot as plt
import traceback


def plot_spectrum(filtered_data, observed_mz, scan_number, window_size, pdf, abundance_threshold=0.0):
    """
    Plots mass spectra for a specific scan and adds the plot to a PDF.
    """
    try:
        # Extract m/z values and relative abundance
        mz_values = filtered_data["m/z Values"].astype(float).values
        relative_abundance = filtered_data["Relative Abundance"].astype(float).values

        # Ensure data is not empty
        if len(mz_values) == 0 or len(relative_abundance) == 0:
            print(f"No valid data for Scan {scan_number}. Skipping plot.")
            return

        # Filter data based on the abundance threshold
        valid_indices = relative_abundance >= abundance_threshold
        filtered_mz = mz_values[valid_indices]
        filtered_abundance = relative_abundance[valid_indices]

        if len(filtered_mz) > 0:
            plt.figure(figsize=(10, 6))
            plt.stem(filtered_mz, filtered_abundance, linefmt='-', markerfmt='o', basefmt=' ')
            plt.xlabel("m/z")
            plt.ylabel("Relative Abundance")
            plt.title(
                f"Scan {scan_number} Spectrum (Filtered {observed_mz - window_size:.2f} - {observed_mz + window_size:.2f} m/z, Abundance ≥ {abundance_threshold})")
            plt.grid(True)
            pdf.savefig()
            plt.close()
            print(f"Plot for Scan {scan_number} added to PDF.")
        else:
            print(
                f"No data points in {observed_mz - window_size:.2f} - {observed_mz + window_size:.2f} m/z range with abundance ≥ {abundance_threshold} for Scan {scan_number}.")

    except Exception as e:
        print(f"Error in plot_spectrum: {e}")
        traceback.print_exc()
