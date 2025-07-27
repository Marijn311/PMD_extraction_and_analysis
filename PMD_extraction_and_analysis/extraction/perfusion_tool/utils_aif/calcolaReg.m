%% ------------------------------------------------------------------------
function [REG] = calcolaReg(y,x,mask)
    % Calculates the irregularity index of the concentration curve for each voxel.
    % The index is calculated by normalizing the area to 1 to avoid penalizing 
    % voxels with high areas.
    % 
    % The formula used to calculate the irregularity index is:
    % CTC = integral((C"(t))^2 dt)
    % where C"(t) is the second derivative of the concentration curve
    %
    % Inputs:
    %   y - 4D concentration data
    %   x - time points
    %   mask - binary mask indicating voxels to process
    
    [nR,nC,nT] = size(y);
    
    % Normalize curves by area under curve (AUC)
    AUC = sum(y,3);
    AUC = AUC + (AUC==1); % Avoid division by zero
    for t=1:nT
        y(:,:,t) = y(:,:,t)./AUC;
    end
    
    % Calculate second derivative
    derivative2 = zeros(nR,nC,nT);
    for r=1:nR
        for c=1:nC
            if mask(r,c)==1
                for t=1:nT
                    if (t>1)&&(t<nT)
                        % Standard case - central difference
                        derivative2(r,c,t) = ((y(r,c,t+1)-y(r,c,t))/(x(t+1)-x(t))-(y(r,c,t)-y(r,c,t-1))/(x(t)-x(t-1)))/(x(t)-x(t-1));
                    elseif t==1
                        % First point - forward difference
                        derivative2(r,c,t) = (y(r,c,t+1)-y(r,c,t))/((x(t+1)-x(t))^2);
                    else
                        % Last point - backward difference
                        derivative2(r,c,t) = (y(r,c,t)-y(r,c,t-1))/((x(t)-x(t-1))^2);
                    end
                end
            end
        end
    end
    
    % Calculate irregularity index by integrating squared second derivative
    derivative2 = derivative2.^2;
    REG = trapz(x,derivative2,3);
    end % calcolaReg
        
       
       