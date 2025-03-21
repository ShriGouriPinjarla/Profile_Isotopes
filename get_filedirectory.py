import os
import json
import mimetypes

def get_excel_and_mzxml_paths(config_path):

    try:
        # Load configuration file
        with open(config_path, 'r') as config_file:
            config = json.load(config_file)

        # Get the base directory containing folders
        base_directory = config["Isotopic_Distribution"]["file_path"][0]
        base_directory = base_directory.strip('"')  # Remove any extra quotes

        if not os.path.exists(base_directory):
            raise ValueError(f"The directory {base_directory} does not exist.")

        # Dictionary to store paths
        folder_file_paths = {}

        # Iterate over each folder in the base directory
        for folder_name in os.listdir(base_directory):
            folder_path = os.path.join(base_directory, folder_name)

            if os.path.isdir(folder_path):
                # Initialize placeholders for Excel and mzXML file paths
                excel_file_path = None
                mzxml_file_path = None

                # Check all files in the folder
                for file_name in os.listdir(folder_path):
                    file_path = os.path.join(folder_path, file_name)
                    mime_type, _ = mimetypes.guess_type(file_path)

                    # Check if the file is an Excel file
                    if mime_type == "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet":
                        excel_file_path = file_path

                    # Check if the file is an mzXML file (specific MIME types may vary)
                    if file_name.endswith('.mzXML') or file_name.endswith('.mzxml'):
                        mzxml_file_path = file_path

                # Add paths only if both files are found
                if excel_file_path and mzxml_file_path:
                    folder_file_paths[folder_name] = (excel_file_path, mzxml_file_path)
                else:
                    # Print warnings if files are missing
                    if not excel_file_path:
                        print(f"No Excel file found in folder: {folder_path}")
                    if not mzxml_file_path:
                        print(f"No mzXML file found in folder: {folder_path}")
            else:

                print(f"Skipping non-directory item: {folder_path}")

        return folder_file_paths

    except Exception as e:
        print(f"Error accessing directories: {e}")
        return {}