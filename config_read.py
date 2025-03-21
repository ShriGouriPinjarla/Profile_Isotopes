import json
from formula_read import *

def validate_modification(modification_str):
    try:
        mass, formula = modification_str.split(':')
        # Check if mass starts with '+' or '-' and is a valid float
        if not (mass.startswith('+') or mass.startswith('-')):
            raise ValueError("Mass must start with '+' or '-'.")
        float(mass)
        parse_formula(formula)
        return True
    except (ValueError, IndexError) as e:
        print(f"Invalid modification format: {modification_str}. Error: {e}")
        return False
def get_config(config_path):
    try:
        with open(config_path, 'r') as file:
            config = json.load(file)

        # Access raw file paths
        file_path = config['Isotopic_Distribution']['file_path']

        # Extract the list of modifications
        modification_objects = config["Isotopic_Distribution"].get("Modification", [])

        # Check if there are any modifications to add
        if not modification_objects:
            print("No modifications to add. Exiting.")
            exit()

        # Load the existing modifications from modifications.py
        with open('modifications.py', 'r') as mod_file:
            mod_lines = mod_file.readlines()

        # Find the line where the modifications dictionary is defined
        db_start = None
        db_end = None
        for i, line in enumerate(mod_lines):
            if 'modifications = {' in line:
                db_start = i
            if db_start is not None and '}' in line:
                db_end = i
                break

        if db_start is not None and db_end is not None:
            # Extract existing modifications to check for duplicates
            existing_modifications = {}
            for i, line in enumerate(mod_lines[db_start + 1:db_end], db_start + 1):
                if ':' in line:
                    key = line.split(':')[0].strip().strip(" '\"")
                    existing_modifications[key] = i


            # Add or replace each new modification in the dictionary
            for modification_str in modification_objects:
                if not validate_modification(modification_str):
                    print(f"Skipping invalid modification: {modification_str}")
                    continue

                # Split the modification string into mass and formula
                mass, formula = modification_str.split(':')

                # Parse the formula into a dictionary
                formula_dict = parse_formula(formula)

                # Create the new modification entry
                new_entry = f"    '{mass}': {formula_dict},\n"

                pattern = re.compile(rf"(\s*'{re.escape(mass)}'\s*:\s*\{{[^}}]*\}},?)")

                modified = False
                for i, line in enumerate(mod_lines):
                    if pattern.search(line):
                        mod_lines[i] = new_entry
                        modified = True
                        print(f"Modification {mass} already exists. Replacing with the new one.")
                        break

                if not modified:
                    # Add the new modification before the closing '}'
                    mod_lines.insert(db_end, new_entry)
                    print(f"Added new modification: {mass}")


            # Write the updated content back to modifications.py
            with open('modifications.py', 'w') as mod_file:
                mod_file.writelines(mod_lines)
            print("Modifications updated successfully!")
        else:
            print("Could not find modifications dictionary in modifications.py.")

        return file_path

    except Exception as e:
        print(f'Error occurred while loading config file: {e}')
        return None
