import pandas as pd
import os
import re
import openpyxl


def extract_scan_number(input_string):

    match = re.search(r'scan=(\d+)', input_string)
    if match:
        return int(match.group(1))
    return None


def extract_sheet1(file_path):

    try:
        # Load the Excel file
        excel_data = pd.ExcelFile(file_path, engine="openpyxl")
        sheet1_data = excel_data.parse(">300score" or "Sheet1")

        # Define required columns and their new names
        required_columns = {
            'Query #:z': 'Query_Number',
            'z': 'Charge_State',
            'Peptide\n< ProteinMetrics Confidential >': 'Peptide_Sequence',
            'Scan #': 'Scan_Number',
            'Observed\nm/z': 'Observed_mz'
        }

        # Check for missing columns
        missing_columns = [col for col in required_columns if col not in sheet1_data.columns]
        if missing_columns:
            raise ValueError(f"Missing columns in 'Sheet1': {missing_columns}")

        # Extract and rename required columns
        extracted_data = sheet1_data[list(required_columns.keys())].copy()
        extracted_data.rename(columns=required_columns, inplace=True)

        # Extract numeric scan number using the extract_scan_number function
        extracted_data['Scan_Number'] = extracted_data['Scan_Number'].astype(str).apply(extract_scan_number)

        # Check if any scan numbers are missing
        if extracted_data['Scan_Number'].isnull().any():
            raise ValueError("Some scan numbers could not be extracted from the 'Scan #' column.")

        # Define output file path
        output_dir = os.path.dirname(file_path)
        output_file = os.path.join(output_dir, "extracted_sheet1_data.csv")

        # Save to CSV format
        extracted_data.to_csv(output_file, index=False)
        print(f"Extracted data saved to: {output_file}")

        return output_file

    except Exception as e:
        print(f"Error during extraction: {e}")
        return None

