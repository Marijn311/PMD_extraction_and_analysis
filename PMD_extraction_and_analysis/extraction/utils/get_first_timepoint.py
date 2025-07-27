import nibabel as nib
import sys

"""
This script is used to extract the first timepoint from a perfusion image.
Effectively taking the first 3D image from a 4D image.
This is necessary to skull strip the perfusion image.
The skull stripping algorithm we use (BET) does not work on 4D images.
Hence the first volume is skull stripped, and this mask is used for all volumes, since making a new mask for every volume is way to slow.
"""

def get_first_timepoint(paths):

    # Extract the perf path (assumed to be the first passed path)
    perf_path = paths[0]  

    # Load the perfusion image	
    img = nib.load(perf_path)
    data = img.get_fdata()
    
    # Extract the first timepoint
    data = data[:,:,:,0]  # Assuming the time dimension is the last one, we take the first timepoint
    img.header.set_data_shape(data.shape) # Update the img header to reflect the new shape
    
    # Save the firdt volume as a new NIfTI file
    first_timepoint_img = nib.Nifti1Image(data, img.affine, img.header)
    new_perf_path = perf_path.replace('_reoriented.nii.gz', '_first_timepoint.nii.gz')
    nib.save(first_timepoint_img, new_perf_path)

    # Update the perf_path in paths
    paths[0] = new_perf_path

    # By printing we can return to the powershell
    print(paths) 
    
if __name__ == "__main__":
    
    # Get paths from the powershell script
    paths = sys.argv[1:]
    
    get_first_timepoint(paths)