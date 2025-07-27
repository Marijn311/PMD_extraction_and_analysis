
% This is the core processing function for Dynamic Susceptibility Contrast (DSC) MRI analysis

% The function handles:
% 1) Input Parameter Validation
% 2) Brain Mask Generation/Validation
% 3) Contrast Concentration Calculation
% 4) Arterial Input Function (AIF) Extraction
% 5) Perfusion Map Calculations


function [cbv,cbf,mtt,cbv_lc,ttp,tmax,mask,aif,conc,s0,fwhm,K2_map]=DSC_mri_core(volumes,te,tr,options,aif,mask,atlas,perf_path)
%add the utils subfolder to the path
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')));

%----------------------------------------------------------------------------------------------------
%----- PART 1: Input Parameter Validation ----------------------
%----------------------------------------------------------------------------------------------------


% Check if volumes is a 4D matrix
if not(ndims(volumes)==4)
    error('Input data "volumes" must be a 4D matrix of a brain DSC-MRI.')
end

% Check if volumes is a double, if not convert it
if not(isa(volumes,'double'))
    volumes=double(volumes);
end

% Check if te and tr are numbers
if not(isnumeric(te)) || not(isnumeric(tr))
    error('Echo time (te) and Repetition time (tr) must be numbers.')
end

% Check if options-struct is provided, otherwise load default
if isempty(options)
    % Load default options. This is also where (display) settings etc are set.
    options=DSC_mri_getOptions; 
else
    % Use provided options
    options=options;
end

% If an external mask is provided, check if it is a double, else convert it
if ~isempty(mask) && ~isa(mask, 'double')
    mask = double(mask);
end

% If an external mask is provided, check if it is a 3D matrix
if ~isempty(mask) && not(ndims(mask)==3)
    error('Input mask must be a 3D matrix.')
end

% If an external mask is provided, check if it is a binary or has values between 0 and 1
if ~isempty(mask) && not(all(all(all(mask==0 | mask==1))))
    error('Input mask must be a binary matrix.')
end

% Check if the atlas is provided, if not give an error
if isempty(atlas)
    error('An atlas is required for rCBF determination.')
end

% make sure the atlas is a double
if ~isa(atlas, 'double')
    atlas = double(atlas);
end


% Print welcome message
if options.display > 0
    disp(' ');
    disp(' ');
    disp(' ');
    disp('  _____   ____    _____         _    _   ____    ___ ')
    disp(' |  _  \ /  __|  /  ___|       | \  / | |  _  \ |   |')
    disp(' | | | | | |__   | |      ___  |  \/  | | |_| |  | | ')
    disp(' | | | | \___ \  | |     |___| | |\/| | |     /  | | ')
    disp(' | |_| |  ___| | | |___        | |  | | | |\ \   | | ')
    disp(' |_____/ |____/  \_____|       |_|  |_| |_| \_\ |___|')
    disp(' ')
    disp('                       TOOLBOX                       ')
    disp('          by Denis Peruzzo & Marco Castellaro        ')
    disp('          Adapted by Marijn Borghouts                ')
    disp(' ')
end

if options.display > 0
    disp('Checking data...');
end

%save the perf_path in the options struct
options.perf_path=perf_path;

% Extract image dimensions and save them in the options struct
volumes = double(volumes);
[nR,nC,nS,nT]=size(volumes);
options.nR=nR;
options.nC=nC;
options.nS=nS;
options.nT=nT;

% Save the acquisition parameters in the options struct
options.te=te;
options.tr=tr;
options.time=0:tr:(nT-1)*tr;

% Display data and acquisition parameters
if options.display > 0
    disp('                      DATA SIZE ');
    disp(['                     Rows - ' num2str(options.nR)]);
    disp(['                  Columns - ' num2str(options.nC)]);
    disp(['                   Slices - ' num2str(options.nS)]);
    disp(['                  Samples - ' num2str(options.nT)]);
    disp(['                Echo time - ' num2str(options.te) ' s']);
    disp(['          Repetition time - ' num2str(options.tr) ' s']);
end


% If an external AIF is provided as input, check if it has the correct dimensions to match the provided data
if ~isempty(aif)
    if not(nT==size(aif.conc,1))
        if not(nT==size(aif.conc,2))
            if size(aif.conc,1)==0 && size(aif.conc,2)==0
                options.aif=1;
            else
                error('AIF dimensions are not compatible to samples data, please provide a compatible AIF.')
            end
        end
    end
end



%----------------------------------------------------------------------------------------------------
%----- PART 2: Brain Mask Generation/Validation -----------------------
%----------------------------------------------------------------------------------------------------

if isempty(mask)
    if not(options.conc) % if no external mask is provided and the input volume is raw singal, generate a mask automatically
        [mask]=DSC_mri_mask(volumes,options);
    else % if no external mask is provided and the input volume is contrast concentration, throw an error
        error('If provided input volume contains concentrations instead of raw signal, you need to also provide a brain mask.')
    end

elseif size(mask,1)==0 %if a mask seems to be provided but is empty, generate a new mask
    if not(options.conc)
        [mask]=DSC_mri_mask(volumes,options);
    end

else %if the mask is provided and not empty, check if it has the correct dimensions to match the provided volume data
    if options.aif.enable
        
        if not(nR==size(mask,1)) || not(nC==size(mask,2)) || not(nS==size(mask,3))
            disp('WARNING! Mask dimensions do not match data dimensions.')
            disp('Mask has dimensions:')
            disp(size(mask))
            disp('Data has dimensions:')
            disp([nR,nC,nS])
            
            error('Mask dimensions do not match data dimensions.')
            % disp('Rescaling mask to data dimensions...')
            % temp = mask;
            % clear mask;
            % mask = imresize3(temp, [nR, nC, nS], 'nearest');
            % mask = mask > 0.5;
        end
        [mask]=DSC_mri_mask_only_aif(volumes,mask,options);
    end

end

%----------------------------------------------------------------------------------------------------
%----- Part 2.5: Detecting / Correcting Patient motion  ------------------------------
%----------------------------------------------------------------------------------------------------

[volumes] = DSC_mri_patient_motion(volumes,mask,options);

%----------------------------------------------------------------------------------------------------
%----- Part 3: Contrast Concentration Calculation ------------------------------
%----------------------------------------------------------------------------------------------------
if options.conc % if the input volume is contrast concentration, assign volume variable to conc variable
    conc=volumes;
    options.aif.nSlice=aif.sliceAIF;
else
    [conc,s0,bolus]=DSC_mri_conc(volumes,mask.data,options);
    conc = real(conc);
end


%----------------------------------------------------------------------------------------------------
%----- Part 4: Arterial Input Function (AIF) extraction ------------------------------
%----------------------------------------------------------------------------------------------------
if options.aif.enable % if in options the AIF calculation is enabled, call the AIF extraction function
    % Ensure mask has an aif field
    if ~isfield(mask, 'aif')
        mask.aif = mask.data;
    end
    [aif]=DSC_mri_aif(conc,mask.aif,options, atlas);
end


%----------------------------------------------------------------------------------------------------
%----- Part 5: Perfusion Map Calculations ------------------------------
%----------------------------------------------------------------------------------------------------

if options.display > 0
    disp(' ')
    disp('Calculating maps...')
end
    
% Call function to calculate CBV map
[cbv]=DSC_mri_cbv(conc,aif.fit.gv,mask,options, atlas);

% Call function to calculate CBV leackage corrected map #this one does not make the value relative by reference.
[cbv_lc,~,K2_map,~,~]=DSC_mri_cbv_lc(conc,aif.fit.gv,mask,bolus,options);

% Call function to calculate CBF map
[cbf]=DSC_mri_cbf(conc,aif.fit.gv,mask,options, atlas, s0, volumes);

% Call function to calculate MTT map
[mtt]=DSC_mri_mtt(cbv_lc,cbf,options);

% Call function to calculate TTP map and FWHM
[ttp]=DSC_mri_ttp(conc,mask.data,atlas,options);
[fwhm]=DSC_mri_fwhm(conc,mask.data,options);

% Call function to calculate Tmax map
[tmax]=DSC_mri_tmax(conc,aif.fit.gv,mask,atlas,options, s0, volumes);

% Output refinements 
mask=mask.data;

end