% DSC_mri_toolbox 
% This is a free and open-source software (FOSS) for the processsing of 
% Dynamic Susceptibility Contrast MRI (DSC-MRI) data into perfusion maps.

%------ Load the dataset to be analyzed ---------------------------------

function DSC_main_demo(perf_path, atlas_path, mask_path)
    % Remove leading and trailing spaces from the paths
    perf_path = strtrim(perf_path);
    atlas_path = strtrim(atlas_path);
    mask_path = strtrim(mask_path);

    % If paths have more than one \ or /, remove them
    perf_path = regexprep(perf_path, '[\\/]+', '/');
    atlas_path = regexprep(atlas_path, '[\\/]+', '/');
    mask_path = regexprep(mask_path, '[\\/]+', '/');

    % Load the images
    DSC_info   = niftiinfo(perf_path);
    DSC_volume = niftiread(DSC_info);

    atlas_info = niftiinfo(atlas_path);
    atlas = niftiread(atlas_info);

    mask_info = niftiinfo(mask_path);
    mask = niftiread(mask_info);

    TE = 0.03; % in seconds
    TR = 1.41;  % in seconds


    % ------ Perform quantification ------------------------------------------ 
    % INPUTS:
    %   volumes - 4D matrix of DSC-MRI data [nR x nC x nS x nT] : nR = rows, nC = columns, nS = slices, nT = time points
    %   te      - Echo time (s)
    %   tr      - Repetition time (s) 
    %   options - Structure containing processing options (optional)
    %   aif     - Arterial input function (optional, default is automatic extraction)
    %   mask    - Binary mask of brain tissue (optional, default is automatic generation)
    %   atlas   - Atlas used to guide automatic AIF extraction

    % OUTPUTS:
    %   cbv     - 3D matrix with rCBV values
    %   cbf     - Struct with 3D matrices of rCBF values for each deconvolution method
    %   mtt     - Struct with 3D matrices of MTT values for each deconvolution method
    %   cbv_lc  - 3D matrix with leackage corrected rCBV values.
    %   ttp     - 3D matrix with leackage corrected Time to Peak values
    %   mask    - 3D matrix with the generated brain mask
    %   aif     - Struct with AIF extracted with clustering algorithm
    %   conc    - 4D matrix with pseudo-concentration values
    %   s0      - 3D matrix with S0 estimates from pre-contrast images
    %   fwhm    - Full Width at Half Maximum 
    %   K2_map  - Estimated K2 leakage map 

    % If the options struct is not provided, the toolbox will take the default
    if ~exist('options', 'var')
        options = [];
    end

    % If no AIF is provided, the toolbox will extract it automatically
    if ~exist('aif', 'var')
        aif = [];
    end

    % If no mask is provided, the toolbox will generate one
    if ~exist('mask', 'var')
        mask = [];
    end

    % if no atlas is provided, the toolbox will give an error at a later point
    if ~exist('atlas', 'var')
        atlas = [];
    end

    % -------Call the main function with the given inputs --------------------
    [cbv,cbf,mtt,cbv_lc,ttp,tmax,mask,aif,conc,s0,fwhm,K2_map]=DSC_mri_core(DSC_volume,TE,TR, options, aif, mask, atlas, perf_path);

    % % ------ View Results ----------------------------------------------------
    % DSC_mri_show_results(cbv_lc,cbf,mtt,ttp,tmax,mask,aif,conc,s0);

    % ------ Save all 3D spatial results in NIFTI format ---------------------
    nifti_template = DSC_info;
    nifti_template.Datatype='double';
    nifti_template.BitsPerPixel=64;
    nifti_template.ImageSize=size(cbv);
    nifti_template.PixelDimensions=nifti_template.PixelDimensions(1:3);


    methods = fieldnames(cbf);
    for k = 1:length(methods)
        deconv_method = methods{k};
        deconv_folder = fullfile(fileparts(perf_path), deconv_method);
        mkdir(deconv_folder)
     
        niftiwrite(cbv,fullfile(deconv_folder, 'CBV.nii'),nifti_template, 'Compressed', true)
        niftiwrite(cbv_lc,fullfile(deconv_folder, 'CBV_LC.nii'),nifti_template, 'Compressed', true)
        niftiwrite(ttp,fullfile(deconv_folder, 'TTP.nii'),nifti_template, 'Compressed', true)
        niftiwrite(mask,fullfile(deconv_folder, 'mask.nii'),nifti_template, 'Compressed', true)
        niftiwrite(s0,fullfile(deconv_folder, 'S0.nii'),nifti_template, 'Compressed', true)

        niftiwrite(cbf.(deconv_method).map,fullfile(deconv_folder, 'CBF.nii'),nifti_template, 'Compressed', true)
        niftiwrite(mtt.(deconv_method),fullfile(deconv_folder, 'MTT.nii'),nifti_template, 'Compressed', true)
        niftiwrite(tmax.(deconv_method),fullfile(deconv_folder, 'Tmax.nii'),nifti_template, 'Compressed', true)
    end

end
