% This function loads the default options for the DSC_mri_core function
% and returns them in a structure.

function [optionsOUT]=DSC_mri_getOptions()


%----------------------------------------------------------------------------------------------------
% DISPLAY OPTIONS
%----------------------------------------------------------------------------------------------------
optionsOUT.display=0; % 2:off, 1:notify (text), 2:notify (images), 3:debug (be carefull, this will open a lot of figures, if there are many slices in the volume)
optionsOUT.waitbar=0; % 0:off, 1:on

optionsOUT.strel_disk = 2; %This is the radius of the disk used for erosion of the reference region mask	

warning('off', 'all');
 
% DATA PREPARATION OPTIONS
optionsOUT.mask.npixel=300; 
% represents the minimum number of pixels of a connected component that is
% used as a threshold for automatic mask generation

optionsOUT.conc=0; 
% 0: the provided data in "volumes" is raw signal, 1: the provided data is contrast agent concentrations

optionsOUT.S0.nSamplesMin=5; 
% Minimum number of initial scans on which to calculate S0

optionsOUT.S0.nSamplesMax=12; 
% Maximum number of initial scans on which to calculate S0

optionsOUT.S0.thresh=0.05; 
% Threshold used to choose the instant of tracer appearance

optionsOUT.min_mask_size=800;

%----------------------------------------------------------------------------------------------------
% OPTIONS FOR AIF IDENTIFICATION PHASE
%----------------------------------------------------------------------------------------------------
optionsOUT.aif.enable= 1; 
% 0: does not calculate the AIF, 1: calculates the AIF

optionsOUT.aif.ricircolo= 1;
% 0: does not account for recirculation, 1: fits the recirculation

optionsOUT.aif.nSlice= 0;
% Slice on which to search for the AIF (0: lets the operator select the slice)

optionsOUT.aif.semiasseMaggiore= 0.3500; 
% Size of the major axis for the search area

optionsOUT.aif.semiasseMinore= 0.1500; 
% Size of the minor axis for the search area

optionsOUT.aif.pArea= 0.4000; 
% Percentage of voxels discarded due to AUC

optionsOUT.aif.pTTP= 0.4000; 
% Percentage of voxels discarded due to TTP

optionsOUT.aif.pReg= 0.0500; 
% Percentage of voxels discarded due to the regularity of the trend

optionsOUT.aif.diffPicco= 0.0400; 
% Threshold to decide whether to select the cluster based on the peak or TTP

optionsOUT.aif.nVoxelMax= 10; 
% maximum number of voxels chosen for the AIF

optionsOUT.aif.nVoxelMin= 5; 
% minimum number of voxels chosen for the AIF

% Correction of the formula for calculating concentration from the signal
% in the case of AIF calculation
optionsOUT.qr.enable= 0; % 0: does not apply the correction, 1: applies the correction
optionsOUT.qr.b= 5.7400e-004;
optionsOUT.qr.a= 0.0076;
optionsOUT.qr.r= 0.0440;



%----------------------------------------------------------------------------------------------------
% OPTIONS FOR DECONVOLUTION METHODS
%----------------------------------------------------------------------------------------------------
optionsOUT.deconv.SVD.threshold= 0.2; %0.2 in original % SVD threshold. This threshold determine how much of the singular values are kept.  it helps in regularizing the inversion of the matrix to avoid amplifying noise.
optionsOUT.deconv.SVD.residual= 1; % 0: does not save residues, 1: saves residues

optionsOUT.deconv.cSVD.threshold= 0.1; % cSVD threshold (0.1 is for data obtained at 1.5T)
optionsOUT.deconv.cSVD.residual= 1; % 0: does not save residues, 1: saves residues

optionsOUT.deconv.oSVD.OIthres = 0.035;    % 10% threshold as in Ostergaard and Calamante
optionsOUT.deconv.oSVD.OIcounter = 1;
optionsOUT.deconv.oSVD.residual= 1; % 0: does not save residues, 1: saves residues


%----------------------------------------------------------------------------------------------------
% TO ADD PARAMETERS FOR STABLE SPLINE
%----------------------------------------------------------------------------------------------------
optionsOUT.deconv.SS.residual = 1;

% Only use one method for deconvolution
optionsOUT.deconv.method={'SVD';}; % Methods to apply for perfusion calculation
% optionsOUT.deconv.method={'SVD'; 'cSVD'; };
% available option {'SVD';'cSVD';'oSVD'} (SVD doesn’t consider the AIF delay effect. oSVD is time consuming. )

%----------------------------------------------------------------------------------------------------
% PROPORTIONALITY CONSTANTS
%----------------------------------------------------------------------------------------------------
optionsOUT.par.kh= 1;
optionsOUT.par.rho= 1;
optionsOUT.par.kvoi= 1;

