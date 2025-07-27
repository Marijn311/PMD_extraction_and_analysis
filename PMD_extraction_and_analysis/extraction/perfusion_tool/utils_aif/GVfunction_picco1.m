


    %% ------------------------------------------------------------------------
    function [GV] = GVfunction_picco1(p,options)
        % Calculates the gamma variate function for the first/main peak
        % The gamma variate function is defined by:
        %
        % GV(t) = A*((t-t0)^alpha)*exp(-(t-t0)/beta)
        %
        % Parameters: p=[t0 alpha beta A]
        % t0 - time offset
        % alpha - shape parameter
        % beta - scale parameter
        % A - amplitude
        
        t0 = p(1);    % Time offset
        alpha = p(2);  % Shape parameter
        beta = p(3);   % Scale parameter
        A = p(4);      % Amplitude
        
        % Initialize output vector
        nT = length(options.time);
        GV = zeros(nT,1);
        
        % Calculate gamma variate function at each timepoint
        for i=1:nT
            t = options.time(i);
            % Only calculate for times after t0 to avoid complex numbers
            if t>t0
                GV(i) = A*((t-t0)^alpha)*exp(-(t-t0)/beta);
            end
        end
        end % GVfunction_picco1
        
            
            
            