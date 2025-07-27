


    %% ------------------------------------------------------------------------
    function [recirculation] = GVfunction_ricircolo(p,options)
        % Calculates the gamma variate function describing the recirculation phase
        % The complete model is:
        %
        % FP(t) = A*((t-t0)^alpha)*exp(-(t-t0)/beta)
        % recirculation(t) = FP(t-td) convolved with K*exp(-t/tao)
        %
        % Parameters are passed in order:
        % p = [t0 alpha beta A td K tao]
        %   t0 - time offset
        %   alpha - shape parameter
        %   beta - scale parameter
        %   A - amplitude
        %   td - time delay
        %   K - scaling factor
        %   tao - exponential decay constant
        %
        % Note: Since the formula involves convolution, calculations are done on a finer
        % time grid than the final output grid for better accuracy
        
        % Extract parameters
        if length(p) ~= 7
            disp('Breakpoint: Parameter vector p must have length 7');
            keyboard; % This will pause execution and allow you to inspect variables
        end
        
        t0 = p(1);    % Time offset
        alpha = p(2);  % Shape parameter
        beta = p(3);   % Scale parameter
        A = p(4);      % Amplitude
        td = p(5);     % Time delay
        K = p(6);      % Scale factor
        tao = p(7);    % Decay constant
        
        % 1) Define fine time grid for convolution calculations
        TR = options.time(2)-options.time(1); % Original time step
        Tmax = max(options.time);
        Tmin = min(options.time);
        nT = length(options.time);
        
        TRfine = TR/10; % Use 10x finer grid
        tGrid = Tmin:TRfine:2*Tmax;
        nTfine = length(tGrid);
        
        % 2) Calculate functions needed for recirculation
        % FP(t) = A*((t-t0)^alpha)*exp(-(t-t0)/beta)
        % disp(t) = exp(-t/tao) 
        % recirculation(t) = K * [FP(t-td) convolved with disp(t)]
        
        % Initialize vectors
        peak2 = zeros(nTfine,1);   % Delayed first pass curve
        dispersion = zeros(nTfine,1); % Exponential dispersion
        
        for i=1:nTfine
            t = tGrid(i);
            
            % Calculate delayed first pass curve
            if t>t0+td
                peak2(i) = K*((t-t0-td)^alpha)*exp(-(t-t0-td)/beta);
            end
            
            % Calculate exponential dispersion
            dispersion(i) = exp(-t/tao);
        end
        
        % 3) Convolve components to get recirculation on fine grid
        recirculation_fine = TRfine.*filter(peak2,1,dispersion);
        
        % 4) Interpolate to get values at requested timepoints
        recirculation = interp1(tGrid,recirculation_fine,options.time);
        end % GVfunction_ricircolo
        
        