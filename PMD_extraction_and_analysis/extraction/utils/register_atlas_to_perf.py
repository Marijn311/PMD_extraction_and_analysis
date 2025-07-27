import ants
import sys
import os

"""
This script is used to register the atlas and template to the perfusion image.
It uses the ANTsPy library to perform the registration.
It uses the SyN method for non-linear registration. 
Since the atlas and template are not the same brain as the patient so we want to allow some local deformations.
"""

def register_atlas_to_perf(paths):
    
    perf_path = paths[0] # Asummed to be the first path
    atlas_path = paths[-2]  # Assumed to be the second last path
    template_path = paths[-1]   # Assumed to be the last path

    # Load the images
    atlas = ants.image_read(atlas_path)
    template = ants.image_read(template_path)
    perf = ants.image_read(perf_path)
    
    # Register the template to the perfusion image. 
    # This time we use the SyN method since the patient and the atlas are not from the same patient.
    # Hence, we want to allow non-linear transformations.
    registration_object = ants.registration(fixed=perf, moving=template, type_of_transform='SyN')

    # Apply the forward transformation to both the template and the atlas, to get them in the perf space
    registered_template = ants.apply_transforms(fixed=perf, moving=template, transformlist=registration_object['fwdtransforms'], interpolator='nearestNeighbor')
    registered_atlas = ants.apply_transforms(fixed=perf, moving=atlas, transformlist=registration_object['fwdtransforms'], interpolator='nearestNeighbor')

    # Save the registered atlas
    registered_atlas_path = os.path.join(os.path.dirname(perf_path), "atlas_registered.nii.gz")
    registered_template_path = os.path.join(os.path.dirname(perf_path), "template_registered.nii.gz")

    ants.image_write(registered_atlas, registered_atlas_path)
    ants.image_write(registered_template, registered_template_path)

    # Update the paths to the registered images
    paths[-2] = registered_atlas_path 
    paths[-1] = registered_template_path
    
    # By printing, we return to the powershel script 
    print(paths)

if __name__ == "__main__":
    
    # Get paths from powershell script
    paths = sys.argv[1:]
    
    register_atlas_to_perf(paths)
