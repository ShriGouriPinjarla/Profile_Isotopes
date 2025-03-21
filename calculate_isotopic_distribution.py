from collections import Counter
from collections import defaultdict
from isotopes_db import *
import re
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal, getcontext
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')  # Use 'Agg' for non-interactive plots
import matplotlib.pyplot as plt

# Define proton mass
PROTON_MASS = Decimal('1.007276')  # Mass of a proton in Da

# Define adduct masses
ADDUCT_MASSES = {
    'H': Decimal('1.00728'),  # Proton
    'Na': Decimal('22.989769'),  # Sodium
    'K': Decimal('38.963708')  # Potassium
}

def gaussian_broadening(masses, abundances, resolution=0.01, sigma=0.0005, target_max_height=100.0):

    # Convert Decimal masses to float
    masses = [float(mass) for mass in masses]

    # Create a mass grid
    min_mass = min(masses) - 0.5
    max_mass = max(masses) + 0.5
    mass_grid = np.arange(min_mass, max_mass, resolution)
    broadened_abundances = np.zeros_like(mass_grid)

    # Apply Gaussian broadening
    for mass, abundance in zip(masses, abundances):
        gaussian = norm.pdf(mass_grid, loc=mass, scale=sigma)
        broadened_abundances += abundance * gaussian

    # Calculate the dynamic scaling factor
    current_max_height = np.max(broadened_abundances)
    scaling_factor = target_max_height / current_max_height

    # Apply the scaling factor
    broadened_abundances = broadened_abundances * scaling_factor

    return mass_grid, broadened_abundances

def convolve_distributions_profile(dist1, dist2, threshold=Decimal('1e-5'), resolution=0.01):

    masses1, abundances1 = dist1
    masses2, abundances2 = dist2

    convolved_masses = []
    convolved_abundances = []

    for i in range(len(masses1)):
        for j in range(len(masses2)):
            convolved_masses.append(float(masses1[i]) + float(masses2[j]))
            convolved_abundances.append(float(abundances1[i]) * float(abundances2[j]))

    # Generate a mass grid based on resolution
    min_mass = min(convolved_masses)
    max_mass = max(convolved_masses)
    mass_grid = np.arange(min_mass, max_mass + resolution, resolution)

    # Aggregate abundances into the mass grid
    aggregated_abundances = np.zeros_like(mass_grid)

    for mass, abundance in zip(convolved_masses, convolved_abundances):
        index = int((mass - min_mass) / resolution)
        aggregated_abundances[index] += abundance

    # Prune the distribution to keep only significant peaks
    max_abundance = np.max(aggregated_abundances)
    significant_indices = aggregated_abundances >= max_abundance * float(threshold)

    pruned_masses = mass_grid[significant_indices]
    pruned_abundances = aggregated_abundances[significant_indices]

    # Normalize abundances to ensure they sum to 1
    total_abundance = np.sum(pruned_abundances)
    normalized_abundances = pruned_abundances / total_abundance

    return pruned_masses, normalized_abundances

def isotopic_distribution_theoretical(formula, charge=1, adduct='H', resolution=0.00001, sigma=0.01):

    masses = [Decimal('0.0')]
    abundances = [Decimal('1.0')]
    current_distribution = (masses, abundances)

    for element, count in formula.items():
        element_masses = elements[element]['masses']
        element_abundances = elements[element]['abundances']
        element_distribution = (element_masses, element_abundances)

        # Convolution logic
        binary_representation = bin(count)[2:]
        intermediate_distribution = element_distribution
        result_distribution = current_distribution

        for bit in reversed(binary_representation):
            if bit == '1':
                result_distribution = convolve_distributions_profile(
                    result_distribution, intermediate_distribution, resolution=resolution
                )
            intermediate_distribution = convolve_distributions_profile(
                intermediate_distribution, intermediate_distribution, resolution=resolution
            )

        current_distribution = result_distribution

    # Apply Gaussian broadening
    masses, abundances = map(list, current_distribution)
    broadened_masses, broadened_abundances = gaussian_broadening(
        masses, abundances, resolution=0.001, sigma=0.0005
    )

    # Adjust masses for adduct and charge state
    adduct_mass = ADDUCT_MASSES.get(adduct, PROTON_MASS)  # Default to proton if adduct is unknown
    adjusted_masses = [
        Decimal(str(mass)) / charge
        for mass in broadened_masses
    ]

    # Convert Decimal masses to float for plotting
    adjusted_masses = [float(mass) for mass in adjusted_masses]

    return adjusted_masses, broadened_abundances

def plot_graph(broadened_masses, broadened_abundances, formula, charge):

    # Plot the isotopic distribution
    plt.figure(figsize=(12, 6))

    # Smooth curve
    plt.plot(broadened_masses, broadened_abundances, color="green",
             label=f"{formula} (Charge {charge}+)", linewidth=1)

    # Labels and title
    plt.xlabel("m/z")
    plt.ylabel("Relative Abundance (%)")
    plt.title("Isotopic Distribution")
    plt.grid(visible=True, linestyle="--", alpha=0.5)

    # Save the plot to a file
    plt.savefig("isotopic_distribution.png", dpi=300, bbox_inches="tight")

    # Extract m/z range (min and max of broadened masses)
    min_mz = min(broadened_masses) - 5
    max_mz = max(broadened_masses) + 5

    return (min_mz, max_mz)