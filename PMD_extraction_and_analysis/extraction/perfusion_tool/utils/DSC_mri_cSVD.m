function [res_csvd]=DSC_mri_cSVD(conc,aif,mask,options)
% last modification: Denis Peruzzo 07/06/2010

% Function from DSC_mri package - DSC_mri_cbf
% Author: Castellaro Marco - University of Padova - DEI
%
% Calculates parametric maps of Cerebral Blood Flow (CBF) for a subject
% using SINGULAR VALUE DECOMPOSITION BLOCK-CIRCULANT method with truncation
% for deconvolution
%
% Input parameters:  - conc (4D Matrix), contains the DSC concentration
%                     time curves for all voxels.
%                   - aif, concentration time curve at the site
%                     chosen as arterial
%                   - mask (3D Matrix), contains the matrix used
%                     to mask the brain volume to analyze
% Options is the struct that contains method options, the significant
% ones are:
%
% options.deconv.svd.threshold - truncation threshold percentage relative to
%                               the maximum eigenvalue, in Wu et al.
%                               Calamante is set to 10%
%
% options.deconv.SVD.residual - if 1 also outputs the 4D matrix
%                              of residuals, otherwise they are not
%                              calculated

% CBF: The CBF for each voxel is determined as the maximum  value of the residue function vettRes.
% Tmax: Tmax for each voxel is determined as the time at which the residue function vettRes reaches its  maximum value.

if options.display > 1
    disp('    Method: cSVD');
    disp(['    Threshold: ' num2str(100*options.deconv.cSVD.threshold,'%3d') '%']);
end

% 1) Create matrix G
nTpad=2*options.nT;
columnG=zeros(nTpad,1);
columnG(1)=aif(1);
columnG(options.nT)=(aif(options.nT-1)+4*aif(options.nT))/6;
columnG(options.nT+1)=aif(options.nT)/6;
for k=2:(options.nT-1)
    columnG(k)=(aif(k-1)+4*aif(k)+aif(k+1))/6;
end
rowG=zeros(1,nTpad);
rowG(1)=columnG(1);
for k=2:nTpad
    rowG(k)=columnG(nTpad+2-k);
end

G=toeplitz(columnG,rowG);

% 2) Apply SVD to calculate inverse G
[U,S,V]=svd(G);

eigenV=diag(S);
threshold=options.deconv.cSVD.threshold*max(eigenV);    % threshold of 10% as in Ostergaard and Calamante
newEigen=zeros(size(eigenV));
for k=1:length(eigenV);
    if eigenV(k)>=threshold;
        newEigen(k)=1/eigenV(k);
    end
end

Ginv=V*diag(newEigen)*(U');

% 3) Apply Ginv to calculate the residue function and CBF of each
%    voxel.

% Initialize the output structure since the code below only fills the pixels that are in the mask
res_csvd.map=zeros(options.nR,options.nC,options.nS);
res_csvd.tmax=zeros(options.nR,options.nC,options.nS); % Initialize tmax field


if options.deconv.SVD.residual
    res_csvd.residual=zeros(options.nR,options.nC,options.nS,nTpad);
end

if options.waitbar
    hw_csvd=waitbar(0,'Computing CBF by cSVD');
end
for r=1:options.nR
    if options.waitbar
        waitbar((r-1)/(options.nR),hw_csvd);
    end
    for c=1:options.nC
        for s=1:options.nS
            if mask(r,c,s)
                % Calculate the residue function
                vettConc=zeros(nTpad,1);
                vettConc(1:options.nT)=reshape(conc(r,c,s,:),options.nT,1);
                vettRes=(1/options.tr)*Ginv*vettConc;
                

                % Due to a combination of patient motion, and oscillations due to the ill-posed nature of deconv.
                % the residue function may have multiple peaks. To avoid artifacts, we assume tha the true peak 
                % is in the first 10 time points.
                
                % Take the highest peak in the first 10 timepoints
                [maxVal, maxIdx] = max(vettRes(1:10));
               
                % the CBF for each voxel is determined as the maximum  value of the residue function vettRes.
                res_csvd.map(r,c,s)=maxVal;
                % Tmax for each voxel is determined as the time at which the residue function vettRes reaches its  maximum value.
                res_csvd.tmax(r,c,s) = maxIdx; 

                if options.deconv.cSVD.residual
                    res_csvd.residual(r,c,s,:)=vettRes;
                end
            end
        end
    end
end

if options.waitbar
    delete(hw_csvd)
end

end