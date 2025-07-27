
    
%% ------------------------------------------------------------------------
function [GVparameters, cv_est_parGV] = fitGV_picco2(data,weights,first_peak_params,options_DSC)
    % Fits the recirculation (second) peak using a gamma variate function
    % The complete model is:
    %
    % FP(t) = A*((t-t0)^alpha)*exp(-(t-t0)/beta)  [First peak]
    % c(t) = FP(t) + FP(t-td) conv K*exp(-t/tao)  [First peak + Recirculation]
    %
    % Parameters: p=[td K tao]
    % td - time delay for recirculation
    % K - scaling factor for recirculation
    % tao - exponential decay constant
    %
    % The fitting is done on the residuals after subtracting the first peak fit
    
    % Set up optimization options
    optimOptions = optimset('lsqnonlin');
    optimOptions.Display = 'none';
    optimOptions.MaxFunEvals = 10000;
    optimOptions.MaxIter = 10000;
    optimOptions.TolX = 1e-8;
    optimOptions.TolFun = 1e-8;
    optimOptions.LargeScale = 'on';
    
    % Ensure all inputs are column vectors
    if size(options_DSC.time,1)==1
        options_DSC.time = options_DSC.time';
    end
    if size(data,1)==1
        data = data';
    end
    if size(weights,1)==1
        weights = weights';
    end
    if size(first_peak_params,1)==1
        first_peak_params = first_peak_params';
    end
    
    % Calculate residuals by subtracting first peak fit
    first_peak = GVfunction_picco1(first_peak_params,options_DSC);
    recirculation_data = data - first_peak;
    
    % Set uniform weights for recirculation fit
    recirculation_weights = ones(length(recirculation_data),1);
    
    % Reduce weight of data before main peak to avoid fitting artifacts
    % Find where concentration drops below 40% of peak after the maximum
    weight_cutoff = min([find(data>0.4*max(data),1,'last'), ...
                        3+find(data==max(data))]);
    recirculation_weights(1:weight_cutoff) = 1;
    
    % Initialize parameter estimates
    % Only consider data after concentration drops below 40% of peak
    init_data = recirculation_data;
    init_data(1:weight_cutoff) = 0;
    init_data(init_data<0) = 0;
    
    % Find time delay (td) as difference between first peak t0 and 
    % recirculation onset (10% of max)
    [recirculation_peak,recirculation_peak_time] = max(init_data);
    recirculation_start = find(init_data(1:recirculation_peak_time)< ...
                                (0.1*max(init_data)),1,'last');
                            
                                
    td_init = options_DSC.time(recirculation_start)-first_peak_params(1);
    
    % Handle edge case where no recirculation is found
    % Sometimes init_data only contains zeros. I guess this happens when there is no clear recirculation peak.
    if isempty(recirculation_start)
        disp('Warning: No recirculation peak found for this slice. so we skip this slice');
        td_init = 0;
    end
    
    
    % Initial estimate for exponential decay constant
    tao_init = 40;
    
    % Scale factor K chosen to match peak heights
    recirculation_init = GVfunction_ricircolo([first_peak_params; td_init; 1; tao_init], ...
                                            options_DSC);
    K_init = max(init_data)./max(recirculation_init);
    
    % Initial parameter vector
    p = [td_init; K_init; tao_init];
    
    % Display initial fit if requested
    if options_DSC.display > 2
        h = figure();
        plot(options_DSC.time,data,'ko',options_DSC.time,first_peak,'k-', ...
                options_DSC.time,GVfunction_ricircolo([first_peak_params; p],options_DSC),'g-')
        title('Recirculation fit - initial values')
        set(h, 'CloseRequestFcn', @my_closefcn);
    end
    
    % Iterative optimization
    cycleFlag = true;
    nCycle = 0;
    while cycleFlag
        % Set parameter bounds
        ub = p.*[10; 10; 10];
        lb = p./[10; 10; 10];
        
        nCycle = nCycle+1;
        [p, resNorm, residuals, exitFlag,OUTPUT,LAMBDA,JACOBIAN] = ...
            lsqnonlin(@objFitGV_picco2, p, lb, ub, optimOptions, recirculation_data, ...
                        recirculation_weights,first_peak_params,options_DSC);
        
        % Exit if converged or max iterations reached
        if (nCycle>=4)||(exitFlag>0)
            cycleFlag = false;
        end
    end
    
    GVparameters = p';
    
    % Calculate parameter uncertainties from covariance matrix
    J = JACOBIAN;
    covp = inv(J'*J);
    var = diag(covp);
    sd = sqrt(var);
    cv_est_parGV = (sd./p*100)';
    
    % Display final fit if requested  
    if options_DSC.display > 2
        figure(h);
        hold on
        plot(options_DSC.time,GVfunction_ricircolo([first_peak_params; p],options_DSC),'r-')
        title('Recirculation final fit')
        set(h, 'CloseRequestFcn', @my_closefcn);
    end
    end % fitGV_picco2
    
    
    
    