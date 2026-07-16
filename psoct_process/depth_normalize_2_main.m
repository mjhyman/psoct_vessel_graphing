%% Main script for Stephan's depth_normalize_2 function
% This script imports a volume, applies the function "depth_normalize_2,"
% and then saves the output to the same directory.

clear; clc; close all;

%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Truncate path to reach top-level directory (psoct_vessel_graphing)
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Initialize data paths
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
% Subfolder containing data
subdir = '/dist_corrected/volume/ref/';
% Filename of OCT volume
volname = 'ref.tif';
% Subject ID list
subid = {'AD_8790'};

%% Import the OCT volume
vol = fullfile(dpath, subid, subdir, volname);
vol = TIFF2MAT(string(vol));

%% Perform normalization with adaptive method
% load reference stack - used to measure mean intensity fluctuation
iref = fullfile(dpath, subid, subdir, 'ref_intensity.tif');
iref = TIFF2MAT(string(iref));
% perform normalization
voln = depth_intensity_normalize(vol,iref);
% Save output
fout = strcat('ref_4ds_norm_alt.tif');
fout = string(fullfile(dpath, subid, subdir, fout));
segmat2tif(voln,fout);

%% Save normalized output as NIFTI
fout = fullfile(dpath, subid, subdir, 'ref_4ds_norm_alt.nii');
save_mri(voln,fout{:},[0.012,0.012,0.012],'float',0);

%% Save mask as NIFTI
% load the mask volume
mask = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/AD_8790/dist_corrected/volume/mask.tif';
mask = TIFF2MAT(mask);
mask = logical(mask);
% save to nifti
fout = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/AD_8790/dist_corrected/volume/mask.nii';
save_mri(mask,fout,[0.012,0.012,0.012],'uchar',0);

%% Perform normalization with Stephan's method
% Set minimum threshold
th = 60;
voln = depth_normalize_2(vol,50,th);
fout = strcat('ref_4ds_norm_',num2str(th),'.tif');
fout = string(fullfile(dpath, subid, subdir, fout));
segmat2tif(voln,fout);
