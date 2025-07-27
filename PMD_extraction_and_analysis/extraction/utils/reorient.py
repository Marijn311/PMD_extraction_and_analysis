import nibabel as nib
import sys

    
"""
This script is used to reorient images to RAS orientation. 
It help when viewing the images that they all have the same orientation.
This also gives the registration algorithms a better starting position.
"""

def reorient(paths):
    """
    Reorient the input by modifying the affine matrix of the NIfTI files.
    We do not modify the data itself so no resampling or interpolation is done.
    """
    paths_to_return = []

    for path in paths:
        
        # Load the image
        img = nib.load(path)

        # Reorient the image to RAS orientation (canonical)
        img_reoriented = nib.as_closest_canonical(img)
      
        # Define the new path for the reoriented image
        new_path = path.replace('.nii.gz', '_reoriented.nii.gz')
     
        # Save the reoriented image
        nib.save(img_reoriented, new_path)

        # Append the new path to the list
        paths_to_return.append(new_path)

    # By printing we return to the powershell script
    print(paths_to_return)


if __name__ == "__main__":
    
    # Get the paths from the powershell script
    paths = sys.argv[1:]
    
    reorient(paths)