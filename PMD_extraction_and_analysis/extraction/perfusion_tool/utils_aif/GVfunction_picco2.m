
    
    
   
%% ------------------------------------------------------------------------
function [GV] = GVfunction_picco2(p,options)
    % Calculates complete gamma variate function with recirculation
    % The model includes first pass and recirculation:
    %
    % FP(t) = A*((t-t0)^alpha)*exp(-(t-t0)/beta)
    % recirculation(t) = FP(t-td) convolved with K*exp(-t/tao)
    %
    % Parameters in order:
    % p = [t0 alpha beta A td K tao]
    %
    % Note: Convolution is done on finer time grid for accuracy
    
    % Calculate first pass curve
    FP = GVfunction_picco1(p(1:4),options);
    
    % Calculate recirculation curve
    recirculation = GVfunction_ricircolo(p,options);
    
    % Combine first pass and recirculation
    GV = FP + recirculation;
    end % GVfunction_picco2
    
       
    