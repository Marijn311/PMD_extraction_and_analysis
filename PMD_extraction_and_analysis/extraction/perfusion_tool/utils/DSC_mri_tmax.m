% C(t) = CBF · AIF(t) ⊗ R(t)
% Tmax is the time to the maximum of the tissue residue function (R).

% Tmax perfusion parameter primarily reflects the bolus delay between the site of the arterial input function (AIF) and the tissue. 

% A drawback of Tmax is that calculation of this perfusion metric requires selection of an AIF (for deconvolution) and that the nature of
% the deconvolution algorithm renders the Tmax perfusion maps very sensitive to even minor changes in the shape of the AIF. 
% So TTP may be more stable, but Tmax considers the aif and can therefore more accurate.



function [tmax]=DSC_mri_tmax(conc,aif,mask,atlas,options, s0, volumes)

    %UPDATE IF YOU INTEND TO ADD A NEW METHOD
    method={'SVD';'cSVD';'oSVD'};

    % Validate requested method exists
    for alg=1:size(options.deconv.method,1)
        if ~ismember(options.deconv.method{alg,:}, method)
            str_error = strjoin(method, ' or ');
            error(['Deconvolution method not recognized. Use: ' str_error])
        end
    end


    for alg=1:size(options.deconv.method,1)
        curr_method = options.deconv.method{alg,:};
        
        % Calculate tmax using specified method 
        % Tmax for each voxel is determined as the time at which the residue function vettRes reaches its abs maximum value.
        % The residue function is the result of the deconvolution.
        % So Tmax is calculated in the Deconvolution methods
        switch curr_method
            case 'SVD'
                residual_function = DSC_mri_SVD(conc,aif,mask.data,options, s0, volumes);
                tmax.svd = residual_function.tmax;
             
            case 'cSVD'
                residual_function = DSC_mri_cSVD(conc,aif,mask.data,options);
                tmax.csvd = residual_function.tmax;
            
            case 'oSVD'
                assert(false, 'i am not using osvd');
        end
    end

    

end


