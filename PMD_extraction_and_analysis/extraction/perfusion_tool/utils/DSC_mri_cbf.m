function [cbf]=DSC_mri_cbf(conc,aif,mask,options, atlas, s0, volumes)
%Calculates the parametric maps of Cerebral Blood Flow (CBF) for a subct
%
%Input parameters - conc (4D Matrix), contains the DSC concentration trends of all voxels.
%                  - aif, trend of concentrations in the site chosen as arterial
%                  - mask (3D Matrix), contains the matrix used to mask the brain volume to be analyzed
%
%Options is the struct that contains the method options, the significant ones are:
%
%The calculation of the CBF parameter requires performing a deconvolution operation, in the toolbox some methodologies suitable for this purpose are provided: SVD and cSVD (Block-circulant version).
%
%options.deconv.method - Must be a cell array containing the names of the algorithms with which you intend to perform the analysis. Eg. SVD cSVD or SS
%
%                       For each method, a struct with the characteristic parameters of the method itself must be inserted in the options, as in this example:
%
%                       options.deconv.<method name>.<par_1>
%                       options.deconv.<method name>.<...>
%                       options.deconv.<method name>.<par_n>
%
%options.time - Represents the time vector of the DSC exam (each sample represents the acquisition of an entire brain volume)
%
%options.display - level 1 Shows the progress of the processing,
%                  level 2 Also gives information on the parameters set for the algorithms used
%
%Output parameters:
%cbf - several sub-structs, one for each method chosen to be used,
%      each sub-struct contains a field map that distinguishes the calculated cbv map, for example:
%      cbf.<method name>.map
%
%      residual, 4D matrix that contains the residuals (it is an optional field that must be requested in the options: parameter options.deconv.<method name>.residual)
%      for example:
%      cbf.<method name>.residual


% It is possible to add a deconvolution method, to do so it is necessary to update the method variable, adding the identifying string and a case for calling the new method, furthermore a function must be written that implements the method and called DSC_mri_<method_name>

%%


% CBF is determined as the absolute max value of the residue function



%UPDATE IF YOU INTEND TO ADD A NEW METHOD
method={'SVD';'cSVD';'oSVD'};

if options.display > 0
    disp('   CBF');
end

% Validate requested method exists
for alg=1:size(options.deconv.method,1)
    if ~ismember(options.deconv.method{alg,:}, method)
        str_error = strjoin(method, ' or ');
        error(['Deconvolution method not recognized. Use: ' str_error])
    end
end

for alg=1:size(options.deconv.method,1)
    curr_method = options.deconv.method{alg,:};
    
    % Calculate CBF using specified method
    switch curr_method
        case 'SVD'
            cbf.svd = DSC_mri_SVD(conc,aif,mask.data,options, s0, volumes);
           
        case 'cSVD'
            cbf.csvd = DSC_mri_cSVD(conc,aif,mask.data,options);
          
        case 'oSVD'
            cbf.osvd = DSC_mri_oSVD(conc,aif,mask.data,options);
           
    end 
end
end
