function [cbv]=DSC_mri_cbv(conc,aif,mask,options, atlas)
% Function of the DSC_mri package - DSC_mri_cbv
% Author: Castellaro Marco - University of Padova - DEI
%
% Calculates the parametric maps of Cerebral Blood Volume for a subject
%
% Input parameters: conc (4D Matrix) containing the DSC concentration
% profiles of all voxels.
% Options is the struct that contains the method options, the significant ones are:
%
% options.time - Represents the time vector of the DSC exam (each
%                sample represents the acquisition of an entire brain volume)
%
% options.par.kh - Parameter representing the dependence on hematocrit
%                  of the Cerebral Blood Volume (CBV), by default set to
%                  one, in this case, relative estimates of the
%                  CBV parameter are obtained.
%
% options.par.rho - Parameter representing the dependence on the density
%                   of blood of the Cerebral Blood Volume (CBV), by default
%                   set to one, in this case, relative estimates
%                   of the CBV parameter are obtained.
%
% Output parameters:
% cbv - (3D matrix) containing the calculated parametric map


% The CBV is proportional to the area under the concentration curve (conc) for each voxel, 
% normalized by the area under the arterial input function curve (aif). 
% The parameters options.par.kh and options.par.rho are used to scale the CBV values.



%CBV is determined as the area under the concentration curve divided by the area under the AIF curve, multiplied by some constants



if options.display > 0
    disp('   CBV');
end

cbv=zeros(options.nR,options.nC,options.nS);

if options.waitbar
    hw_cbv=waitbar(0,'Calculating CBV');
end
for s=1:options.nS
    % CBV is determined as a constant (par/rho). (times the mask to remove the background)
    % times the area under the concentration curve divided by the area under the AIF curve.
    
    % the rho anf kh are constants, These constants are both 1, so 1/1=1 so the constants dont effect the result.


    %the area under the concentration curve is a 256x256 matrix, Which makes sense because then the cbv will be differnet for each voxel.
    %the area under the aif curve is a constant, so the the aif is just used to scale the cbv map, as well.
    
    

    % aif is a 1x80 double, so there is only 1 value per time point. 
    % This is the conctrast concentration in the selected voxels over time?.
    % When the aif is changed the resulting cbv map changes as well. But it all changes with the same factor.
    % So the cbv map doesnt change in shape, just in scale. However when the mean of the cbv was 0, then after scaling the mean is still 0.
    % So I guess the cbv mean can be changed through the conc function.
    

    % % Define the range of pixels to average around the pixel (130,130,s)
    % range = 0:0;
    % avg_signal = zeros(length(options.time), 1);
    % count = 0;
    % 
    % for dx = range
    %     for dy = range
    %         for dz = range
    %             if (130+dx > 0 && 130+dx <= options.nR) && (130+dy > 0 && 130+dy <= options.nC) && (s+dz > 0 && s+dz <= options.nS)
    %                 pixel = conc(130+dx, 130+dy, s+dz, :);
    %                 avg_signal = avg_signal + pixel(:);
    %                 count = count + 1;
    %             end
    %         end
    %     end
    % end
    % 
    % avg_signal = avg_signal / count;

    % % Plot the averaged concentration over time
    % figure
    % plot(options.time, avg_signal)
    % title('Averaged Concentration over time')
    % xlabel('Time')
    % ylabel('Concentration')
    
    factors = (options.par.kh/options.par.rho)*mask.data(:,:,s);
    integral_conc = trapz(options.time,conc(:,:,s,:),4);
    integral_aif = trapz(options.time,aif);
    cbv(:,:,s) = factors.*(integral_conc./integral_aif);
    
    if s == 1 && options.display > 0
        disp("aif integral: " + integral_aif)
    end 
   
    if options.waitbar
        waitbar(s/options.nS,hw_cbv)
    end
end

