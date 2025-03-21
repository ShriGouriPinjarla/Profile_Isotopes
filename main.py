from config_read import *
import os
import json
import glob
from extract_mzxml import *
from get_filedirectory import *
from extract_excel import *
from process_file22222 import *
from process_mzxml import *
import time


def main():
    total_start_time = time.time()
    try:
        config_path = "config.json"

        file_path = get_config(config_path)
        print(f"File Path: {file_path}")

    except Exception as e:
        print(f"No file_path: {e}")

    try:
        # Process each directory in the file path
        for directory in file_path:
            # Get Excel and mzXML file paths for each folder
            folder_paths = get_excel_and_mzxml_paths(config_path)

            # Print the retrieved paths
            for folder, paths in folder_paths.items():
                excel_path, mzxml_path = paths
                print(f"Folder: {folder}")
                print(f"  Excel Path: {excel_path}")
                print(f"  mzXML Path: {mzxml_path}")

                output_excel = extract_sheet1(excel_path)
                if output_excel:
                    print(f"Cleaned Excel data saved at: {output_excel}")

                output_mzxml = process_mzxml(mzxml_path)

                window_size = 10  # Define the m/z window size as needed
                abundance_threshold = 0.15
                process_mzxml_excel(output_mzxml, output_excel, window_size, abundance_threshold, test_mode=True,
                                    test_rows=5)
    except Exception as e:
            print(f"Error in raw file: {e}")

    total_end_time = time.time()

    total_time = total_end_time - total_start_time
    print(f"Total time to run the entire code: {total_time:.4f} seconds")

if __name__ == "__main__":
     main()




