function [res_svd]=DSC_mri_SVD(conc,aif,mask,options, s0, volumes)
% Function of the DSC_mri package - DSC_mri_cbf
% Author: Castellaro Marco - University of Padova - DEI
%
% Calculates the parametric maps of Cerebral Blood Flow (CBF) for a subject
% and uses the SINGULAR VALUE DECOMPOSITION method with truncation for deconvolution
%
% Input parameters - conc (4D Matrix), contains the DSC concentration trends of all voxels.
%                   - aif, concentration trends in the site chosen as arterial
%                   - mask (3D Matrix), contains the matrix used to mask the brain volume to be analyzed
% Options is the struct that contains the method options, the significant ones are:
%
% options.deconv.svd.threshold - truncation threshold percentage referred to
%                                the maximum eigenvalue, in Ostergaard and
%                                Calamante et al. it is fixed at 20%
%
% options.deconv.SVD.residual - if 1 also produces the 4D matrix of residuals as output,
%                               otherwise they are not calculated








% CBF: The CBF for each voxel is determined as the maximum  value of the residue function vettRes.
% Tmax: Tmax for each voxel is determined as the time at which the residue function vettRes reaches its  maximum value.





if options.display > 1
    disp('    Method: SVD');
    disp(['    Threshold: ' num2str(100*options.deconv.SVD.threshold,'%3d') '%']);
end

% 1) Create the matrix G
aifVett=zeros(options.nT,1);
aifVett(1)=aif(1);
aifVett(options.nT)=aif(options.nT);

for k=2:(options.nT-1)
    aifVett(k)=(aif(k-1)+4*aif(k)+aif(k+1))/6;
end

G=toeplitz(aifVett,[aifVett(1) zeros(1,options.nT-1)]);

% 2) Apply SVD to calculate the inverse of G
[U,S,V]=svd(G);

eigenV=diag(S);
threshold=options.deconv.SVD.threshold*max(eigenV);    % 10% threshold in Ostergaard and Calamante
newEigen=zeros(size(eigenV));
for k=1:length(eigenV);
    if eigenV(k)>=threshold;
        newEigen(k)=1/eigenV(k);
    end
end

Ginv=V*diag(newEigen)*(U');

% 3) Apply Ginv to calculate the residual function and the CBF of each voxel.

% Initialize the output structure since the code below only fills the pixels that are in the mask
res_svd.map=zeros(options.nR,options.nC,options.nS);
res_svd.tmax=zeros(options.nR,options.nC,options.nS); % Initialize tmax field

if options.deconv.SVD.residual
    res_svd.residual=zeros(options.nR,options.nC,options.nS,options.nT);
end

if options.waitbar
    hw_svd=waitbar(0,'Computing CBF by SVD');
end

counter = 0;


for r=1:options.nR
    if options.waitbar
        waitbar((r-1)/(options.nR),hw_svd);
    end
    for c=1:options.nC
        for s=1:options.nS
            if mask(r,c,s)
                
                
                
                % Calculate the residual function
                vettConc=reshape(conc(r,c,s,:),options.nT,1);
                vettRes=(1/options.tr)*Ginv*vettConc;
                
                % Due to a combination of patient motion, and oscillations due to the ill-posed nature of deconv.
                % the residue function may have multiple peaks. To avoid artifacts, we assume that the true peak 
                % is in the first 10 time points.
                
                [maxVal, maxIdx] = max(vettRes(1:10));
            
                % the CBF for each voxel is determined as the maximum  value of the residue function vettRes in that voxel. 
                res_svd.map(r,c,s)=maxVal;
                % Tmax for each voxel is determined as the time at which the residue function vettRes reaches its  maximum value in that voxel.
                res_svd.tmax(r,c,s) = maxIdx; 


                % % plot vettRes, 
                % if maxVal > 0.1 && counter < 3
                %     figure; 
                %     hold on;
                %     plot(reshape(conc(r,c,s,:), [], 1));
                %     plot(reshape(volumes(r,c,s,:), [], 1));
                %     title('Concentration vs Time for cases with unusual high CBF');
                %     legend('Concentration', 'Volumes');
                %     hold off;
                %     counter = counter + 1;
                % end
                
                if options.deconv.SVD.residual
                    res_svd.residual(r,c,s,:)=vettRes;
                end
                
            end
            
        end
    end
end

if options.waitbar
    delete(hw_svd)
end


end