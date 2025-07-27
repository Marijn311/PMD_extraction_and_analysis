"""
This is the main script to extract PMDs from dataset with perfusion images.
Here you can set the paths to the dataset, atlas and template.
The script will loop over all folders in the dataset path and perform the following steps:
1. Get the perfusion image and other optional images (such as DWI) from the folder
2. Reorient the images to the same orientation
3. Extract the first timepoint from the perfusion image (used for skull stripping)
4. Skull strip the images using HD-BET
5. Register the atlas to the perfusion image
6. Register the optional images (such as DWI) to the perfusion image
7. Run the perfusion tool from the matlab script, to get perfusion maps from the perfusion image
8. Normalize the images with the healthy white matter
9. Generate the PMD excel files
"""

# Fixed paths
$DATASET_PATH = Join-Path -Path (Resolve-Path "$PSScriptRoot\..\..") "dummy_dataset\dataset"
$ATLAS_PATH = Join-Path -Path (Resolve-Path "$PSScriptRoot\..\..") "dummy_dataset\atlas\HO_atlas.nii.gz"
$TEMPLATE_PATH = Join-Path -Path (Resolve-Path "$PSScriptRoot\..\..") "dummy_dataset\atlas\MNI152lin_T1_1mm_brain.nii.gz"

# Here you can add the filenames of the images you want to process
$IMAGE_FILENAMES = @("perf_nii\perf.nii.gz", "dwi_nii\b1000.nii.gz")

# Loop over all patients in dataset_path: extract the groupname, patient folder name, and the study folder name
$groups = Get-ChildItem -Path $DATASET_PATH -Directory
Write-Host "Found $($groups.Count) groups in $DATASET_PATH"
foreach ($group in $groups) {
    $groupName = $group.Name
    Write-Host "Processing group: $groupName"
    $folders = Get-ChildItem -Path $group.FullName -Directory
    Write-Host "Found $($folders.Count) patient folders in $groupName"
    foreach ($folder in $folders) {
        $folderName = $folder.Name
        Write-Host "Processing patient: $folderName"
        $studyfolder = Get-ChildItem -Path $folder.FullName -Directory
        if ($studyfolder.Count -eq 1) {
            Write-Host "Found studyfolder: $($studyfolder.Name)"
        } else {
            Write-Host "Error: Expected 1 studyfolder, found $($studyfolder.Count)"
            continue
        }
        $studyFolderPath = $studyfolder.FullName

        # Put all the image paths in a list
        # This is the list of image filenames will be the input/output for the python scripts
        $paths = @()
        foreach ($filename in $IMAGE_FILENAMES) {
            $path = Join-Path -Path $studyFolderPath -ChildPath $filename
            $paths += $path
        }
        
        # Reorient all images to the same orientation
        conda activate PMD
        Set-Location utils
        $paths = python reorient.py $paths
        
        # Parse python output to a list again
        $paths = $paths -replace "\[|\]", ""
        $paths = $paths -replace "\\", "/"
        $paths = $paths -replace "'", ""
        $paths = $paths -split ","
        $paths = $paths -replace "^\s+|\s+$", ""
        
        # Save a copy of the reoriented 4D perfusion image. The perfusion toolbox needs it later.
        $PERF_PATH = $paths[0]
        
        Write-Output "Reorientation done"

        # Extract the first timepoint from the perfusion image
        $paths = python get_first_timepoint.py $paths
        
        # Parse python output to a list again
        $paths = $paths -replace "\[|\]", ""
        $paths = $paths -replace "\\", "/"
        $paths = $paths -replace "'", ""
        $paths = $paths -split ","
        $paths = $paths -replace "^\s+|\s+$", ""
        
        Write-Output "First timepoint extracted"

        # Prepare output paths for skull stripping
        $outputPaths = @() 
        foreach ($path in $paths) {
            if ($path -like "*first_timepoint.nii.gz") {
                $outputPath = $path.Replace("first_timepoint.nii.gz", "stripped.nii.gz")
            } else {
                $outputPath = $path.Replace("reoriented.nii.gz", "stripped.nii.gz")
            }
            $outputPaths += $outputPath
        }

        # Parse output paths 
        $outputPaths = $outputPaths -replace "\[|\]", ""
        $outputPaths = $outputPaths -replace "\\", "/"
        $outputPaths = $outputPaths -replace "'", ""
        $outputPaths = $outputPaths -split ","
        $outputPaths = $outputPaths -replace "^\s+|\s+$", ""
        
        # Skull strip the images using HD-BET
        conda activate hdbet
        for ($i = 0; $i -lt $paths.Count; $i++) {
            $inputPath = $paths[$i]
            $outputPath = $outputPaths[$i]
            
            # Run HD-BET on the input path and save the output to the output path
            hd-bet -i $inputPath -o $outputPath -device cpu --disable_tta --save_bet_mask
        }

        # Reset the skullstripped paths as the new paths
        $paths = $outputPaths

        # Add the atlas and template paths to the end of the paths list
        $paths += $ATLAS_PATH
        $paths += $TEMPLATE_PATH
        
        Write-Output "Skull stripping done"

        # Register the atlas to the perfusion image.
        conda activate PMD
        $paths = python register_atlas_to_perf.py $paths

        # Parse python output to a list again
        $paths = $paths -replace "\[|\]", ""
        $paths = $paths -replace "\\", "/"
        $paths = $paths -replace "'", ""
        $paths = $paths -split ","
        $paths = $paths -replace "^\s+|\s+$", ""
        
        Write-Output "Atlas registered to perfusion image"

        # If paths has more than 3 elements (perf, atlas,template)
        if ($paths.Count -gt 3) {
            
            # Register optional images such as dwi to the perfusion image
            $paths = python register_to_perf.py $paths
            
            # Parse python output to a list again
            $paths = $paths -replace "\[|\]", ""
            $paths = $paths -replace "\\", "/"
            $paths = $paths -replace "'", ""
            $paths = $paths -split ","
            $paths = $paths -replace "^\s+|\s+$", ""
            
            Write-Output "Registered alternative images to perf"
        }
        
        # Run the perfusion tool from the matlab script
        Set-Location ..
        Set-Location perfusion_tool
        $mask_path = $paths[0].Replace("stripped.nii.gz", "stripped_bet.nii.gz")    
        $atlas_path = $paths[-2]
        matlab -batch "DSC_main('$PERF_PATH', '$atlas_path', '$mask_path')"
        Set-Location .. 
        Set-Location utils

        Write-Output "Perfusion maps calculated"
        
        # Get the paths to the perfusion maps
        $ParentDirOfPerf = Split-Path -Path $PERF_PATH -Parent
        $SVDPATH = Join-Path -Path $ParentDirOfPerf -ChildPath "svd"
        $CBF_PATH = Join-Path -Path $SVDPATH -ChildPath "CBF.nii.gz"
        $CBV_LC_PATH = Join-Path -Path $SVDPATH -ChildPath "CBV_LC.nii.gz"
        $CBV_PATH = Join-Path -Path $SVDPATH -ChildPath "CBV.nii.gz"
        $MTT_PATH = Join-Path -Path $SVDPATH -ChildPath "MTT.nii.gz"
        $TTP_PATH = Join-Path -Path $SVDPATH -ChildPath "TTP.nii.gz"
        $Tmax_PATH = Join-Path -Path $SVDPATH -ChildPath "Tmax.nii.gz"
        
        # Make a new paths list where the paths are in the order of perf, all the maps from svd folder, all alternatives, atlas, and template
        $new_paths = @($paths[0], $CBF_PATH, $CBV_LC_PATH, $CBV_PATH, $MTT_PATH, $TTP_PATH, $Tmax_PATH)
        for ($i = 0; $i -lt $paths.Count; $i++) {
            $path = $paths[$i]
            if ($i -eq 0) {
                continue
            } else {
                $new_paths += $path
            }
        }
        $paths = $new_paths
        
        # Normalize the perfusion maps
        conda activate PMD
        $paths = python normalize.py $paths

        # Parse python output to a list again
        $paths = $paths -replace "\[|\]", ""
        $paths = $paths -replace "\\", "/"
        $paths = $paths -replace "'", ""
        $paths = $paths -split ","
        $paths = $paths -replace "^\s+|\s+$", ""
        
        Write-Output "Perfusion maps normalized"

        # Get the PMD excel files for this image
        python generate_pmds.py $paths

        Write-Output "PMD excel files generated"

        # Move back to the original directory before processing the next patient
        cd ..
    }
}