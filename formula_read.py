import re
from collections import Counter, defaultdict
from modifications import *
from amino_acids import *
from decimal import Decimal, getcontext
import traceback
import sys


try :
    # Set precision for Decimal operations
    getcontext().prec = 20

    # Define proton mass
    PROTON_MASS = Decimal('1.007276')  # Mass of a proton in Da

    # Define adduct masses
    ADDUCT_MASSES = {
        'H': Decimal('1.007276'),  # Proton
        'Na': Decimal('22.989769'),  # Sodium
        'K': Decimal('38.963708')  # Potassium
    }
    # Water molecule to subtract per peptide bond
    WATER = {'H': 2, 'O': 1}

    def parse_formula(formula_str):

        if not isinstance(formula_str, str):
            raise TypeError(f"Expected a string, but got {type(formula_str).__name__}")

        pattern = r'([A-Z][a-z]*)(\d*)'
        matches = re.findall(pattern, formula_str)  # Find all matches

        formula = {}
        for element, count in matches:
            formula[element] = int(count) if count.isdigit() else 1

        return formula


    def extract_peptide_and_modifications(sequence):

        # Remove flanking residues (e.g., "-.M" and ".F")
        match = re.match(r"^[^A-Za-z]*([A-Za-z].*?)(?:\.[A-Za-z])?$", sequence)
        if match:
            sequence = match.group(1)  # Extract the core sequence

        # Extract modifications (numbers inside square brackets)
        modification_masses = re.findall(r"\[([\+\-]?\d+\.\d+)\]", sequence)

        # Remove modifications from sequence
        sequence = re.sub(r"\[.*?\]", "", sequence)

        return sequence, modification_masses


    def format_modification_key(mod_mass):

        return f"+{float(mod_mass):.3f}" if mod_mass else None

    def calculate_peptide_formula(peptide_sequence, modification_masses):


        element_counts = Counter()

        # Add elements for each amino acid
        for aa in peptide_sequence:
            if aa in amino_acids:
                for element, count in amino_acids[aa].items():
                    element_counts[element] += count
            else:
                raise ValueError(f"Invalid amino acid: {aa}")

        # Subtract water molecules for each peptide bond
        num_peptide_bonds = len(peptide_sequence) - 1  # (n - 1) peptide bonds for n amino acids
        for element, count in WATER.items():
            element_counts[element] -= num_peptide_bonds * count

        # Process modifications and add to the formula
        for mod_mass in modification_masses:
            formatted_mod = format_modification_key(mod_mass)
            if formatted_mod and formatted_mod in modifications:
                for element, count in modifications[formatted_mod].items():
                    element_counts[element] += count
            else:
                print(f"Warning: Modification mass '{mod_mass}' not found in database.")

        # Construct chemical formula as a string
        chemical_formula = "".join(
            f"{element}{count}" for element, count in sorted(element_counts.items()) if count > 0)

        return chemical_formula, dict(element_counts)


except:
    print("Unexpected error:", sys.exc_info()[0])
    print(traceback.print_exc())