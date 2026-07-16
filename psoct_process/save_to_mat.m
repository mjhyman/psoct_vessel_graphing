%% Apply the tissue, white matter, and gray matter masks to segmentation
% This script applies several masks to the segmentation. These are the
% following masks: tissue, white matter (WM), gray matter (GM), sulci, and
% gyri. The names of the masks are: mask_tiss, mask_gm, mask_wm,
% mask_sulci, mask_gyri. The purpose of applying these masks is to analyze
% the various regional characteristics of the brain tissue specimens.
% The masks are all stored in the directory:
% Ann_Mckee_samples_55T/[subID]/dist_corrected/volume/ref/masks/
%
% Then, this script converts the masked components of segmentation into
% graphs before analyzing the graphs.

clear; clc; close all;

%% Add top-level directory of code repository to path
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach parent
% (psoct_human_brain\vasculature\vesSegment)
topdir = mydir(1:idcs(end));
addpath(genpath(topdir));

%% Initialize paralell pool
% Set # threads = # cores for job
NSLOTS = str2num(getenv('NSLOTS'));
maxNumCompThreads(NSLOTS);

%% Import gyri, sulci, wm

% data directory
d = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/AD_8790/dist_corrected/volume/ref/masks/';

% gyri mask
f = fullfile(d,"mask_gyri.tif");
mask = logical(TIFF2MAT(f));
f = fullfile(d,"mask_gyri.mat");
save(f,'mask','-v7.3');

% sulci mask
f = fullfile(d,"mask_sulci.tif");
mask = logical(TIFF2MAT(f));
f = fullfile(d,"mask_sulci.mat");
save(f,'mask','-v7.3');

% WM mask
f = fullfile(d,"mask_wm.tif");
mask = logical(TIFF2MAT(f));
wm_mask = mask;
f = fullfile(d,"mask_wm.mat");
save(f,'mask','-v7.3');

% tissue
f = fullfile(d,'mask_tiss.nii');
mask = MRIread(f,0,0);
mask = logical(mask.vol);
tiss_mask = mask;
f = fullfile(d,"mask_tiss.mat");
save(f,'tiss_mask','-v7.3');

% create GM mask from WM and tissue
mask = logical(wm_mask .* tiss_mask);
f = fullfile(d,"mask_gm.mat");
save(f,'mask','-v7.3');

