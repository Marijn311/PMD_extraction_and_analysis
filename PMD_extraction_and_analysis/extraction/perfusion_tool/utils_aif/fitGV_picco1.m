



%% ------------------------------------------------------------------------
function [GVparameters, cv_est_parGV] = fitGV_picco1(data,weights,options_DSC)
    % Fits first peak of concentration curve with gamma variate function
    %
    % The gamma variate function is:
    % GV(t) = A*((t-t0)^alpha)*exp(-(t-t0)/beta)
    %
    % Parameters: [t0 alpha beta A]
    % t0 - time offset
    % alpha - shape parameter
    % beta - scale parameter 
    % A - amplitude
    %
    % Inputs:
    %   data - concentration time series to fit
    %   weights - fitting weights for each time point
    %   options_DSC - struct containing options
    
    % Set up optimization options
    optimOptions = optimset('lsqnonlin');
    optimOptions.Display = 'none';
    optimOptions.MaxFunEvals = 1000;
    optimOptions.MaxIter = 1000;
    optimOptions.TolX = 1e-4;
    optimOptions.TolFun = 1e-4;
    optimOptions.LargeScale = 'on';
    
    % Initialize parameter estimates
    alpha_init = 5; % Fixed initial shape parameter
    
    % Find t0 as last point before curve exceeds 5% of peak
    [peakValue,peakIndex] = max(data);
    peakTime = options_DSC.time(peakIndex);
    t0_init = options_DSC.time(find(data(1:peakIndex)<=0.05*peakValue, 1, 'last'));
    
    % Estimate beta from time-to-peak relationship
    beta_init = (peakTime-t0_init)./alpha_init;
    
    % Set A to match peak height 
    A_init = peakValue./max(GVfunction_picco1([t0_init; alpha_init; beta_init; 1],options_DSC));
    
    % Set up parameter bounds
    p0 = [t0_init; alpha_init; beta_init; A_init];
    lb = p0.*0.1; % Lower bounds
    ub = p0.*10;  % Upper bounds
    
    % Display initial fit if requested
    if options_DSC.display > 2
        h = figure();
        plot(options_DSC.time,data,'ko',options_DSC.time,GVfunction_picco1(p0,options_DSC),'g-')
        title('First peak fit - initial values')
        set(h, 'CloseRequestFcn', @my_closefcn);
    end
    
    % Ensure inputs are column vectors
    if size(options_DSC.time,1)==1
        options_DSC.time = options_DSC.time';
    end
    if size(data,1)==1
        data = data';
    end
    if size(weights,1)==1
        weights = weights';
    end
    
    % Increase fitting weight near peak
    [~,TTP] = max(data);
    weights(TTP) = weights(TTP)./10;
    weights(TTP-1) = weights(TTP-1)./2;
    
    % Find end of first peak (20% of peak value)
    i = TTP;
    while data(i) > 0.2*data(TTP) 
        if i <= length(data) 
            break %if the threshold of 20 is not reached we take the last point, this can happen in the top slices for example where the aif is definetly wrong. (but this slice wont be selected as the best one anyway)
        end
        i = i+1;
    end
    
    
    
    % Zero data after first peak
    data_peak1 = zeros(size(data));
    data_peak1(1:i) = data(1:i);
    
    weights_peak1 = 0.01+zeros(size(weights));
    weights_peak1(1:i) = weights(1:i);
    
    % Iterative optimization
    cycleFlag = true;
    nCycle = 0;
    p = p0;
    while cycleFlag
        nCycle = nCycle+1;
        [p, resNorm, residuals, exitFlag,OUTPUT,LAMBDA,JACOBIAN] = lsqnonlin(@objFitGV_picco1, p, lb, ub, optimOptions, data_peak1, weights_peak1,options_DSC);
        
        if (nCycle>=4)||(exitFlag>0)
            cycleFlag = false;
        end
    end
    GVparameters = p';
    
    % Calculate parameter uncertainties
    J = JACOBIAN;
    covp = inv(J'*J);
    var = diag(covp);
    sd = sqrt(var);
    cv_est_parGV = (sd./p*100)';
    % Display final fit if requested
    if options_DSC.display > 2
        figure(h);
        hold on
        plot(options_DSC.time,GVfunction_picco1(p,options_DSC),'r-')
        title('First peak final fit')
        set(h, 'CloseRequestFcn', @my_closefcn);
    end
    end % fitGV_picco1
    
    
    
    