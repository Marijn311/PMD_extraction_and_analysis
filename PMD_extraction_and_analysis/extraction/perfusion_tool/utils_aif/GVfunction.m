



%% ------------------------------------------------------------------------
function [GV] = GVfunction(p,options)
    % Calculates complete gamma variate function with recirculation
    % This is an alternative implementation that calculates everything
    % on a fine time grid first, then interpolates to final timepoints
    %
    % The model is:
    % FP(t) = A*((t-t0)^alpha)*exp(-(t-t0)/beta)
    % GV(t) = FP(t) + [FP(t-td) convolved with K*exp(-t/tao)]
    %
    % Parameters:
    % p = [t0 alpha beta A td K tao]
    %   t0 - time offset
    %   alpha - shape parameter
    %   beta - scale parameter
    %   A - amplitude 
    %   td - time delay
    %   K - scale factor
    %   tao - decay constant
    
    % Extract parameters
    t0 = p(1);     
    alpha = p(2);  
    beta = p(3);   
    A = p(4);      
    td = p(5);     
    K = p(6);      
    tao = p(7);    
    
    % 1) Define fine time grid for calculations
    TR = options.time(2)-options.time(1);
    Tmax = max(options.time);
    nT = length(options.time);
    
    TRfine = TR/10;  % 10x finer grid
    tGrid = 0:TRfine:2*Tmax;
    nTfine = length(tGrid);
    
    % 2) Calculate GV components on fine grid
    % Split into main components:
    peak1 = zeros(1,nTfine); % First pass peak
    peak2 = zeros(1,nTfine); % Delayed first pass
    dispersion = zeros(1,nTfine); % Exponential dispersion
    
    for i=1:nTfine
        t = tGrid(i);
        
        % Calculate first pass
        if t>t0
            peak1(i) = A*((t-t0)^alpha)*exp(-(t-t0)/beta);
        end
        
        % Calculate delayed first pass
        if t>t0+td
            peak2(i) = K*((t-t0-td)^alpha)*exp(-(t-t0-td)/beta);
        end
        
        % Calculate dispersion
        dispersion(i) = exp(-t/tao);
    end
    
    % 3) Combine components on fine grid
    recirculation = TRfine.*filter(peak2,1,dispersion);
    conc = peak1 + recirculation;
    
    % 4) Sample at requested timepoints
    GV = zeros(1,nT);
    for i=1:nT
        % Find closest point on fine grid
        [err,pos] = min(abs(tGrid-options.time(i)));
        GV(i) = conc(pos);
        
        % Warn if interpolation error is large
        if err>1
            disp('WARNING: Poor approximation due to grid spacing.')
        end
    end
    end % GVfunction
    