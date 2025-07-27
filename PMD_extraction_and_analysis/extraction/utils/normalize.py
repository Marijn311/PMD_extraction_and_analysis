import nibabel as nib
import numpy as np
from scipy import ndimage
import sys
import os

"""
This script is used to normalize the perfusion maps (CBF/CBV etc) and the alternative images (such as DWI) 
This is done by dividing the data by a reference value.
The reference value is the mean value of a region of healthy white matter.
"""

def normalize(paths):

    perf_path = paths[0]
    atlas_path = paths[-2]
    template_path = paths[-1]

    paths_to_return = []
    paths_to_return.append(perf_path)

    
    # Remove the perfusion path and the atlas and template path from the paths list, so we know which images to normalize.
    # The perfusion path is assumed to be the first path in the list
    # The atlas and template path are assumed to be the last two paths in the list
    paths_to_process = paths[1:-2]

    for path in paths_to_process:
        
        # Load the images
        img_info = nib.load(path)
        img_data = img_info.get_fdata()

        # Check for inf or nan values in the image data
        bad_percentage = (np.isnan(img_data) | np.isinf(img_data)).sum() / img_data.size
        
        # Write to error_logs.txt if more than 1% of the voxels are inf or nan
        if bad_percentage > 0.01:
            if not os.path.exists("error_logs.txt"):
                with open("error_logs.txt", "w") as f:
                    f.write("Error logs: \n")
            with open("error_logs.txt", "a") as f:
                f.write(f"More than 1% of voxels are inf or nan in {path}\n")

        # Set inf and nan values to zero
        img_data[np.isnan(img_data)] = 0
        img_data[np.isinf(img_data)] = 0

        """
        Due to numerical unstability in some of the calculations in the perfusion tool, 
        the perfusion maps can have a few pixels with very high or low values but that are not true inf values.
        These extreme values can greatly affect the mean value of the reference region and thus the normalization process.
        Hence we want to set these extreme values to zero.
        However we don't want to remove "normal" high and low values because that is where the interesting data is.
        We cannot take use the mean or std to remove these values because the mean and std are also affected by the extreme values.
        So I will use the IQR to remove the extreme values.

        The threshold for extreme values is set to x times the IQR.
        """
        
        nonzero_img_data = img_data[img_data != 0]
        percentile_25 = np.percentile(nonzero_img_data, 25)
        percentile_50 = np.percentile(nonzero_img_data, 50)
        percentile_75 = np.percentile(nonzero_img_data, 75)

        low_bound = percentile_25 - (25 * (percentile_50 - percentile_25))
        high_bound = percentile_75 + (25 * (percentile_75 - percentile_50))

        img_data[img_data > high_bound] = 0
        img_data[img_data < low_bound] = 0

        # #plot histogram of img_data
        # from matplotlib import pyplot as plt
        # plt.hist(img_data[img_data!=0].flatten(), bins=100)
        # plt.title(f"Histogram of {path}")
        # plt.show()

        whole_brain_mean = np.mean(img_data[img_data != 0])
        
        # If the mean of the reference area too close to zero, write to error_logs.txt
        if whole_brain_mean < 0.01:
            if not os.path.exists("error_logs.txt"):
                with open("error_logs.txt", "w") as f:
                    f.write("Error logs: \n")
            with open("error_logs.txt", "a") as f:
                f.write(f"The normalization whole-brain-mean in {path} is {whole_brain_mean} \n")

  
        # Normalize the entire image by the mean in the reference region
        normalized_data = img_data / whole_brain_mean # The background will be zero, if zero is divided by a number, it will remain zero.	

        # Save the normalized data to a new NIfTI file
        normalized_img = nib.Nifti1Image(normalized_data, img_info.affine, img_info.header)
        if os.path.basename(path).endswith("registered.nii.gz"):
            new_path = path.replace("registered.nii.gz", "normalized.nii.gz")
        else:
            new_path = path.replace(".nii.gz", "_normalized.nii.gz")
        nib.save(normalized_img, new_path)
        paths_to_return.append(new_path)
    
    # When all paths are processed, also add back the atlas and template paths
    paths_to_return.append(atlas_path)
    paths_to_return.append(template_path)

    # Return to the powershell script by printing
    print(paths_to_return)
    

if __name__ == "__main__":
    
    # Get paths from powershell script
    paths = sys.argv[1:]

    normalize(paths)