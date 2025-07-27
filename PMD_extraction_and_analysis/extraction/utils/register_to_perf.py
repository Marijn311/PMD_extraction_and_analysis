import ants 
import sys

"""
This script is used to register the alternative images (such as a DWI) to the perfusion image.
It uses the ANTsPy library to perform the registration.
It uses the affine method for linear registration.
Since the images are from the same patient, we do not need to allow for local deformations.
"""

def register_to_perf(paths):

    paths_to_return = []
        
    perf_path = paths[0] # Asummed to be the first path
    atlas_path = paths[-2]  # Assumed to be the second last path
    template_path = paths[-1]   # Assumed to be the last path

    perf = ants.image_read(perf_path)

    for path in paths:
    
        # Only keep the paths of the alternative images
        if path == perf_path or path == atlas_path or path == template_path:
            continue

        alt_img = ants.image_read(path)

        # Register the alt_img to the perf image using ANTS. We use an affine transformation since the images are from the same patient. 
        registration_object = ants.registration(fixed=perf, moving=alt_img, type_of_transform='Affine')
        
        # Apply the forward transformation to the alt_img image
        registered_alt = ants.apply_transforms(fixed=perf, moving=alt_img, transformlist=registration_object['fwdtransforms'], interpolator='nearestNeighbor')
        
        # Save the registered alternative image to file
        new_path = path.replace('_stripped.nii.gz', '_registered.nii.gz')
        ants.image_write(registered_alt, new_path)
        
        # Append the new path to the list
        paths_to_return.append(new_path)
    
    
    paths_to_return.insert(0, perf_path) # At the perf images path as the start of the list
    paths_to_return.append(atlas_path) # At the atlas images path as the second last of the list
    paths_to_return.append(template_path) # At the template images path as the last of the list
    
    # Return by printing
    print(paths_to_return)	


if __name__ == "__main__":
    
    # Get paths from powershell script
    paths = sys.argv[1:]
    
    register_to_perf(paths)