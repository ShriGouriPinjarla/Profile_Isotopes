import os
import glob
from extract_mzxml import *
def process_mzxml(mzxml_path):
    import os
    try:
        output_data_np = os.path.join(os.path.dirname(mzxml_path), "ms1_scan_data.csv")
        print(f"Extracting MS1 data from {mzxml_path}...")
        scan_data = extract_ms1_scans(mzxml_path)
        print(f"Saving data to {output_data_np}...")
        save_scan_data(scan_data, mzxml_path)
        return output_data_np

    except Exception as e:
        print(f"Error in process_mzxml: {e}")
        traceback.print_exc()
        return None