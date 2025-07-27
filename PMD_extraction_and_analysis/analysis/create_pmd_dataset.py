import pandas as pd
import os
import xarray as xr

"""
This script loops over all patients in the dataset and loads the PMD excel file per patient.
It then creates a 4D xarray with the an axis for the all image types, the cerebral regions, the statistical metrics, and the patients.

This script assume the following structure of the dataset:

dataset
    group_001
        patient_001
            pmds.xlsx
        patient_002
            pmds.xlsx
        patient_003
            ...
    group_002
        patient_001
            pmds.xlsx
        patient_002
            pmds.xlsx
        patient_003
            ...    
    group_003
        ...
"""

DATASET_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "dummy_dataset", "dataset"))
GROUPS = ["group_001", "group_002"]
PMD_FILE_NAME = 'pmds.xlsx'
OUTPUT_NAME = 'pmd_dataset.nc'

# Function to open an excel file with the PMD data and convert it to an xarray
def load_patient_data(patient_folder_path):
    file_path = os.path.join(patient_folder_path, PMD_FILE_NAME)
    
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None
    
    # Create a list to store DataArrays for each sheet
    sheet_arrays = []
    
    # Extract the sheet names from the excel file
    sheet_names = pd.ExcelFile(file_path).sheet_names
      
    # Load each sheet
    for sheet_name in sheet_names:
    
        patient_data = pd.read_excel(file_path, sheet_name=sheet_name, index_col=0)

        # Convert to DataArray
        sheet_array = xr.DataArray(
            patient_data.values,
            dims=["img_type", "region"],
            coords={"img_type": patient_data.index, "region": patient_data.columns}
        )
        
        # Add metric dimension
        sheet_array = sheet_array.expand_dims("metric").assign_coords(metric=[sheet_name])
        sheet_arrays.append(sheet_array)
            
    # Combine all sheets into a single DataArray
    patient_data = xr.concat(sheet_arrays, dim='metric')
    
    return patient_data



if __name__ == "__main__":
    
    for group in GROUPS:
        data_array = None
        dataset_path = os.path.join(DATASET_PATH, group)
        for root, patient_folders, _ in os.walk(dataset_path):
            for patient_folder in patient_folders:
                for subfolder in os.listdir(os.path.join(root, patient_folder)):
                    patient_folder_path = os.path.join(root, patient_folder, subfolder)
                    patient_data = load_patient_data(patient_folder_path)

                    print(f"Processing patient {patient_folder}")

                    #---------------------------------------------------------------------------------------------------
                    if patient_data is not None:
                        patient_data = patient_data.expand_dims("patient").assign_coords(patient=[patient_folder])
                        
                        if data_array is None:
                            data_array = patient_data
                        else:
                            data_array = xr.concat([data_array, patient_data], dim='patient')
                    #---------------------------------------------------------------------------------------------------
            break 
        data_array.to_netcdf(os.path.join(dataset_path, OUTPUT_NAME))
    
print("Finished!")