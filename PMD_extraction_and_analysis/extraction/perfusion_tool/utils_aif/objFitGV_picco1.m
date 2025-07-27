

%% ------------------------------------------------------------------------
function [out] = objFitGV_picco1(p,data,weights,options)
    % Objective function to minimize for fitGV_picco1
    % This calculates the weighted residuals between the gamma variate model 
    % and the measured data for the first/main peak
    %
    % Inputs:
    %   p - Parameter vector [t0, alpha, beta, A]
    %   data - Measured concentration data
    %   weights - Weighting factors for each timepoint
    %   options - Struct containing time vector
    %
    % Output:
    %   out - Vector of weighted residuals
    
    % Calculate predicted values from gamma variate model
    predicted = GVfunction_picco1(p,options);
    
    % Return weighted residuals
    out = (predicted-data)./weights;
    end % objFitGV_picco1
    
    
    