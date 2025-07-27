


%% ------------------------------------------------------------------------
function [out] = objFitGV_picco2(p,data,weights,first_peak_params,options) 
    % Objective function to minimize for the gamma variate fitting of the recirculation peak
    % Calculates weighted residuals between model prediction and measured data
    %
    % Inputs:
    %   p - Parameters for recirculation peak [td, K, tao]
    %   data - Measured concentration data
    %   weights - Fitting weights for each timepoint
    %   first_peak_params - Parameters from first peak fit [t0,alpha,beta,A] 
    %   options - Struct containing time vector
    %
    % Output:
    %   out - Vector of weighted residuals
    
    % Calculate model prediction by combining first peak and recirculation
    predicted = GVfunction_ricircolo([first_peak_params; p],options);
    
    % Return weighted residuals 
    out = (predicted-data)./weights;
    end % objFitGV_picco2
    