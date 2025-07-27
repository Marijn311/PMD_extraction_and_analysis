import os
import pandas as pd
import nibabel as nib
import numpy as np
from diptest import diptest
import sys

"""
This script generates the perfusion map descripters (PMDs) for each patient.
It calculates the mean, standard deviation, median, interquartile range, skewness, kurtosis and Hartigan's dip test for each image in each region of the atlas.
The results are saved in an Excel file.
This script can also processes non-perfusion maps such as DWI.
"""

def load_and_process_image(path):
    """Load and process a single image."""

    img = nib.load(path)
    data = img.get_fdata().astype(np.float64)
    return data, img

def find_contra_region(region):
    """Find the contra region for a given region."""

       
    CONTRA_MAPPING_DICT = {
        1: 49, 2: 50, 3: 51, 4: 52, 5: 53, 6: 54, 7: 55, 8: 56, 9: 57, 10: 58,
        11: 59, 12: 60, 13: 61, 14: 62, 15: 63, 16: 64, 17: 65, 18: 66, 19: 67, 20: 68,
        21: 69, 22: 70, 23: 71, 24: 72, 25: 73, 26: 74, 27: 75, 28: 76, 29: 77, 30: 78,
        31: 79, 32: 80, 33: 81, 34: 82, 35: 83, 36: 84, 37: 85, 38: 86, 39: 87, 40: 88,
        41: 89, 42: 90, 43: 91, 44: 92, 45: 93, 46: 94, 47: 95, 48: 96, 49: 1, 50: 2,
        51: 3, 52: 4, 53: 5, 54: 6, 55: 7, 56: 8, 57: 9, 58: 10, 59: 11, 60: 12,
        61: 13, 62: 14, 63: 15, 64: 16, 65: 17, 66: 18, 67: 19, 68: 20, 69: 21, 70: 22,
        71: 23, 72: 24, 73: 25, 74: 26, 75: 27, 76: 28, 77: 29, 78: 30, 79: 31, 80: 32,
        81: 33, 82: 34, 83: 35, 84: 36, 85: 37, 86: 38, 87: 39, 88: 40, 89: 41, 90: 42,
        91: 43, 92: 44, 93: 45, 94: 46, 95: 47, 96: 48, 97: 108, 99: 110, 100: 111,
        101: 112, 102: 113, 103: 114, 104: 104, 105: 115, 106: 116, 107: 117, 108: 97,
        110: 99, 111: 100, 112: 101, 113: 102, 114: 103, 115: 105, 116: 106, 117: 107
    }

    contra_region = CONTRA_MAPPING_DICT[region]

    # # Print matches to verify that this is correct
    # region_name = MAPPING.loc[region, "name"]
    # contra_region_name = MAPPING.loc[contra_region, "name"]
    # print(f"{region_name} is paired with {contra_region_name}")
    
    return contra_region



def generate_pmds(paths):
    """Generate the PMDs excel file"""
    

    # Load the atlas
    atlas_path = paths[-2]
    atlas_data, _ = load_and_process_image(atlas_path)
    
    # Remove the first and last two paths (this is the perfusion image and the atlas+template)
    paths_to_process = paths[1:-2] 

    # Create dictionary in which we store the img data for every image in the paths_to_process list
    img_data_dict = {}
    img_names = []
    for path in paths_to_process:
        img_name = os.path.basename(path)
        img_name = img_name.replace("_normalized.nii.gz", "") 
        img_names.append(img_name)
        img_data, _ = load_and_process_image(path)
        img_data_dict[img_name] = img_data

    # Check if all images have the same shape
    shapes = [data.shape for data in img_data_dict.values()] + [atlas_data.shape]
    if not all(shape == shapes[0] for shape in shapes):
        raise ValueError(f"The images do not have the same shape. Shapes: {shapes}")

    # Initialize dataframes to store the results
    atlas_regions = np.unique(atlas_data[atlas_data != 0]) # Exclude region 0 (background)
    df_mean = pd.DataFrame(index=img_names, columns=atlas_regions)
    df_std = pd.DataFrame(index=img_names, columns=atlas_regions)
    df_median = pd.DataFrame(index=img_names, columns=atlas_regions)
    df_iqr = pd.DataFrame(index=img_names, columns=atlas_regions)
    df_skew = pd.DataFrame(index=img_names, columns=atlas_regions)
    df_kurt = pd.DataFrame(index=img_names, columns=atlas_regions)
    df_hart = pd.DataFrame(index=img_names, columns=atlas_regions)
    
    # Loop over all atlas regions and calculate the PMDs for every image	
    for atlas_region in atlas_regions:
        
        # Get the mask for the current region
        mask = atlas_data == atlas_region 
        
        for img in img_names:
            data = img_data_dict[img]
            valid_data = data[mask & (data != 0)]
            if len(valid_data) > 0:
                df_mean.at[img, atlas_region] = valid_data.mean()
                df_std.at[img, atlas_region] = valid_data.std()
                df_median.at[img, atlas_region] = np.median(valid_data)
                df_iqr.at[img, atlas_region] = np.percentile(valid_data, 75) - np.percentile(valid_data, 25)
                if valid_data.std() == 0: #to prevent division by zero
                    df_skew.at[img, atlas_region] = np.nan 
                    df_kurt.at[img, atlas_region] = np.nan
                else:
                    df_skew.at[img, atlas_region] = 3 * (valid_data.mean() - np.median(valid_data)) / valid_data.std() # pearson skewness (tells something about the symmetry of the distribution)
                    df_kurt.at[img, atlas_region] = ((valid_data - valid_data.mean())**4).mean() / valid_data.std()**4 # pearson kurtosis (tells something about the shape of tails of the distribution, long or short tails)
                df_hart.at[img, atlas_region], _ = diptest(valid_data)  # hartigan's dip test (tells something about the deviation from unimodality)
            else: # This might happen when a (small) region does not exist anymore after (a bad) registration
                df_mean.at[img, atlas_region] = np.nan
                df_std.at[img, atlas_region] = np.nan
                df_median.at[img, atlas_region] = np.nan
                df_iqr.at[img, atlas_region] = np.nan
                df_skew.at[img, atlas_region] = np.nan
                df_kurt.at[img, atlas_region] = np.nan
                df_hart.at[img, atlas_region] = np.nan
                

    # Save both every metric to a separate sheet in the same excel file
    pmds_path = os.path.join(os.path.dirname(os.path.dirname(atlas_path)), "pmds.xlsx") #We get the dirname twice to go back one more level in the directory structure	
    with pd.ExcelWriter(pmds_path, engine='openpyxl') as writer:
        df_mean.to_excel(writer, sheet_name="mean")
        df_std.to_excel(writer, sheet_name="std")
        df_median.to_excel(writer, sheet_name="median")
        df_iqr.to_excel(writer, sheet_name="iqr")
        df_skew.to_excel(writer, sheet_name="skew")
        df_kurt.to_excel(writer, sheet_name="kurt")
        df_hart.to_excel(writer, sheet_name="hart")
    
    return pmds_path
   

def append_asymmetry_pmds(pmds_path, paths):

    # Load the pmds excel and make a list which contains all the sheet names (metrics such as mean, std, median, etc.)
    xls = pd.ExcelFile(pmds_path)
    all_metrics = xls.sheet_names

    for metric in all_metrics:
        
        # Load the data from the current metric
        df = pd.read_excel(pmds_path, sheet_name=metric)

        # Set the first column (the ones with the names of the image types) as index
        df = df.set_index(df.columns[0])
        
        # Make a df which is copy of the original df. Than clear all the values except the indexes and column header. This will be used to store the absolute differences between the left and right side of the brain
        df_diff = df.copy()
        df_diff.iloc[:,:] = np.nan
    
        # Loop over all regions. Per region we find the contra lateral region and 
        # calculate the absolute difference between the corresponding regions for each image type
        asymmetry_list = []
        for region in df.columns:

            contra_region = find_contra_region(int(region))
            
            for image_type in df.index:

                # Get the value of cell with region / image_type that we are looking at now and its contra lateral region
                region_value = df.loc[image_type, region]
                contra_region_value = df.loc[image_type, contra_region]

                # Determine the abs difference between the region and its contra lateral region
                abs_diff = abs(region_value - contra_region_value)

                # Save the absolute difference in a list
                asymmetry_list.append((region, contra_region, abs_diff))
            
                # Save the absolute difference in the df_diff dataframe
                df_diff.loc[image_type, region] = abs_diff
                df_diff.loc[image_type, contra_region] = abs_diff
                

        # Open the Excel file in append mode and write df_diff to a new sheet without "Unnamed: 0"
        with pd.ExcelWriter(pmds_path, mode="a", engine="openpyxl", if_sheet_exists="replace") as writer:
            
            # Write the DataFrame to Excel
            df_diff.to_excel(writer, sheet_name=f"{metric}_diff", index=True)
            
            # Clear cell A1 in the header row, this says unnamed: 0
            worksheet = writer.sheets[f"{metric}_diff"]
            worksheet.cell(row=1, column=1).value = None


# Loop over all images (This assumes a specific folder structure)
if __name__ == "__main__":
    
    # Get paths from powershell script
    paths = sys.argv[1:]

    pmds_path = generate_pmds(paths)
    append_asymmetry_pmds(pmds_path, paths)
    

                   