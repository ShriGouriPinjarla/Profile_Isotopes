import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import traceback
from calculate_isotopic_distribution import *
from formula_read import *
from calculate_isotopic_distribution import *
from modifications import *
from amino_acids import *
from isotopes_db import *
from calculate_dense_plot import *

def process_mzxml_excel(output_mzxml, output_path, window_size, abundance_threshold, test_mode=False, test_rows=1):
    try:

        # Validate file paths
        if not os.path.exists(output_path):
            raise FileNotFoundError(f"The file {output_path} does not exist.")
        if not os.path.exists(output_mzxml):
            raise FileNotFoundError(f"The file {output_mzxml} does not exist.")

        # Load the saved CSV data and convert to DataFrames
        spectra_data = pd.read_csv(output_path)
        mzxml_data = pd.read_csv(output_mzxml)

        # Validate required columns
        required_columns = ['Query_Number', 'Charge_State', 'Peptide_Sequence', 'Observed_mz']
        if not all(column in spectra_data.columns for column in required_columns):
            raise ValueError(f"The file {output_path} is missing required columns.")

        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(output_path)
        os.makedirs(output_dir, exist_ok=True)

        # Create a PDF file to save the plots
        output_pdf = os.path.join(output_dir, "interactive_spectra_plots.pdf")
        print(f"Saving plots to: {output_pdf}")

        # Limit the number of rows to process in test mode
        if test_mode:
            spectra_data = spectra_data.head(test_rows)
            print(f"Running in test mode. Processing only the first {test_rows} row(s).")

        with PdfPages(output_pdf) as pdf:
            for _, row in spectra_data.iterrows():
                query_number = row['Query_Number']
                charge_state = row['Charge_State']
                peptide_sequence = row['Peptide_Sequence']
                observed_mz = row['Observed_mz']

                try:
                    adduct = 'H'
                    # Extract the valid peptide sequence and modifications
                    peptide, modifications = extract_peptide_and_modifications(peptide_sequence)

                    # Compute the chemical formula with modifications
                    formula, element_counts = calculate_peptide_formula(peptide, modifications)

                    # Print the result
                    print("Peptide Sequence:", peptide_sequence)
                    print("Chemical Formula:", formula)

                    masses_t, abundances_t = isotopic_distribution_theoretical(
                        element_counts, charge=charge_state, adduct=adduct, resolution=0.0001, sigma=0.002
                    )

                    # Generate a smooth curve for the isotopic distribution
                    broadened_masses, broadened_abundances = gaussian_broadening(
                        masses_t, abundances_t, resolution=0.001, sigma=0.002
                    )


                except:
                    print("Error in theoretical file:", sys.exc_info()[0])
                    print(traceback.print_exc())

                    continue

                # Calculate m/z range
                mz_min = int(np.floor(observed_mz - window_size))
                mz_max = int(np.ceil(observed_mz + window_size))


                # Initialize variables to store the best scan data
                best_scan_number = None
                best_filtered_mz = None
                best_filtered_abundances = None
                max_total_abundance = 0

                zoom_mz = None
                zoom_abundances = None

                # Iterate through all scans in the mzXML data
                for scan_number, scan_data in mzxml_data.groupby('Scan_Number'):
                    try:
                        # Extract m/z and relative abundance data manually
                        mz_values_str = scan_data['m/z Values'].values[0]
                        abundance_values_str = scan_data['Relative Abundances'].values[0]

                        if not isinstance(mz_values_str, str) or not isinstance(abundance_values_str, str):
                            print(f"Warning: Invalid m/z values for Scan {scan_number}")
                            continue

                        # Preprocess the m/z Values column
                        mz_values_list = [float(x) for x in mz_values_str.strip('[]').split(',') if x.strip()]
                        mz_array = np.array(mz_values_list)

                        # Preprocess the Relative Abundances column
                        abundance_values_list = [float(x) for x in abundance_values_str.strip('[]').split(',') if x.strip()]
                        abundance_array = np.array(abundance_values_list)

                        # Check if arrays are empty after conversion
                        if mz_array.size == 0 or abundance_array.size == 0:
                            print(f"Warning: No valid m/z values found in range for Scan {scan_number}")
                            continue

                        # Filter data within m/z range and above abundance threshold
                        valid_indices = (mz_array >= mz_min) & (mz_array <= mz_max) & (abundance_array >= abundance_threshold)
                        filtered_mz = mz_array[valid_indices]
                        filtered_abundances = abundance_array[valid_indices]

                        # Calculate total abundance for the filtered data points
                        total_abundance = np.sum(filtered_abundances)

                        # Update the best scan data if the current scan has higher total abundance
                        if total_abundance > max_total_abundance:
                            max_total_abundance = total_abundance
                            best_scan_number = scan_number
                            best_filtered_mz = filtered_mz
                            best_filtered_abundances = filtered_abundances

                            # Find the densest region
                            zoom_mz_min, zoom_mz_max = find_densest_region(best_filtered_mz, best_filtered_abundances)

                            # Find the nearest data points to zoom_mz_min and zoom_mz_max
                            zoom_mz_min = best_filtered_mz[np.argmin(np.abs(best_filtered_mz - zoom_mz_min))]
                            zoom_mz_max = best_filtered_mz[np.argmin(np.abs(best_filtered_mz - zoom_mz_max))]

                            zoom_indices = (best_filtered_mz >= zoom_mz_min) & (best_filtered_mz <= zoom_mz_max)
                            zoom_mz = best_filtered_mz[zoom_indices]
                            zoom_abundances = best_filtered_abundances[zoom_indices]

                    except Exception as e:
                        print(f"Error processing Scan {scan_number}: {e}")
                        traceback.print_exc()

                if best_scan_number is not None:
                    # Create a figure with 2x2 subplots
                    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

                    # Set the title for the entire page
                    fig.suptitle(f"Peptide: {peptide_sequence}", fontsize=16)

                    # Plot theoretical distribution
                    ax1.stem(broadened_masses, broadened_abundances, linefmt='b-', markerfmt=' ', basefmt=' ',
                             label='Theoretical')
                    ax1.set_xlabel("m/z")
                    ax1.set_ylabel("Relative Abundance")
                    ax1.set_title("Theoretical Distribution")
                    ax1.grid(True)
                    ax1.legend()

                    ax2.stem(zoom_mz, zoom_abundances, linefmt='r-', markerfmt=' ', basefmt=' ', label='Experimental')

                    # Add labels, title, and grid
                    ax2.set_xlabel("m/z")
                    ax2.set_ylabel("Relative Abundance")
                    ax2.set_title("Experimental Distribution in Densest Region")
                    ax2.grid(True)
                    ax2.legend()

                    # Plot overlapping distributions
                    ax3.stem(broadened_masses, broadened_abundances, linefmt='b-', markerfmt=' ', basefmt=' ',
                             label='Theoretical')
                    ax3.stem(best_filtered_mz, best_filtered_abundances, linefmt='r-', markerfmt=' ', basefmt=' ',
                             label='Experimental')
                    ax3.set_xlabel("m/z")
                    ax3.set_ylabel("Relative Abundance")
                    ax3.set_title("Overlapping Distributions")
                    ax3.grid(True)
                    ax3.legend()

                    # Calculate similarity scores
                    cosine_similarity_score = calculate_similarity_score(broadened_masses, broadened_abundances,
                                                                         best_filtered_mz, best_filtered_abundances)
                    peak_pattern_similarity_score = calculate_peak_pattern_similarity(broadened_masses,
                                                                                      broadened_abundances,
                                                                                      best_filtered_mz,
                                                                                      best_filtered_abundances)
                    rms_difference = calculate_rms_difference(broadened_masses, broadened_abundances, best_filtered_mz,
                                                              best_filtered_abundances)

                    # Display similarity scores and RMS difference
                    ax4.text(0.1, 0.8, f"Cosine Similarity: {cosine_similarity_score:.2f}", fontsize=12)
                    ax4.text(0.1, 0.6, f"Peak Pattern Similarity: {peak_pattern_similarity_score:.2f}", fontsize=12)
                    ax4.text(0.1, 0.4, f"RMS Difference: {rms_difference:.2f}", fontsize=12)
                    ax4.axis('off')  # Hide axes for the text box

                    # Save figure to the interactive PDF
                    pdf.savefig(fig)
                    plt.close(fig)

                    print(f"Plot for Scan {best_scan_number} added to PDF.")
                else:
                    print(
                        f"No data points in {observed_mz - window_size:.2f} - {observed_mz + window_size:.2f} m/z range "
                        f"with abundance ≥ {abundance_threshold} for Peptide {peptide_sequence}."
                    )
        # Verify PDF creation
        if os.path.exists(output_pdf):
            print(f"PDF successfully saved: {output_pdf}")
        else:
            print("Error: PDF was not saved!")

    except Exception as e:
        print(f"Error in process_mzxml_excel: {e}")
        traceback.print_exc()