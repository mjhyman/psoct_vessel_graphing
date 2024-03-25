%% Convert mask.mat and ref.mat to NIFTI format
% Each subdirectory ".../[subjectID]/dist_corrected/volume/ref contains the
% c-scans and mask c-scans in the format ref.mat and mask.mat. The mask
% must be corrected in freeview to remove false positives/negatives. This
% script will iterate through these files for each subject and create a
% NIFTI for each.

clear; clc; close all;

%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders

% Start in current directory
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

%% Initialize data path for linux or personal machine (debugging)

%%% Local machine
if ispc
    dpath = ['C:\Users\mack\Documents\BU\Boas_Lab\psoct_data_and_figures'...
            '\test_data\Ann_Mckee_samples_10T\'];
%%% Computing cluster (SCC)
elseif isunix
    % Path to top-level directory
    dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';
    % Set # threads = # cores for job
%     NSLOTS = str2num(getenv('NSLOTS'));
%     maxNumCompThreads(NSLOTS);
% 
%     % Check to see if we already have a parpool, if not create one with
%     % our desired parameters
%     poolobj = gcp('nocreate');
%     if isempty(poolobj)
% 	    myCluster=parcluster('local');
% 	    % Ensure multiple parpool jobs don't overwrite other's temp files
% 	    myCluster.JobStorageLocation = getenv('TMPDIR');
% 	    poolobj = parpool(myCluster, NSLOTS);
%     end
end

% Subfolder containing data
subdir = '/dist_corrected/volume/ref';
% Filename to parse (this will be the same for each subject)
fbase = 'ref';
%%% Complete subject ID list for Ann_Mckee_samples_10T
subid = {'AD_10382','AD_20832', 'AD_20969',...
         'AD_21354', 'AD_21424',...
         'CTE_6489','CTE_6912',...
         'CTE_7019','CTE_8572','CTE_7126',...
         'NC_6839','NC_6974','NC_8653',...
         'NC_21499','NC_301181'};
subid = {'AD_10382'};

for ii = 1:length(subid)
% parfor (ii = 1:length(subid), NSLOTS)
    %% Load "maskec.nii" from the subfolder
    % Filenames
    maskname = 'maskec.nii';
    % Define entire filepath 
    maskpath = fullfile(dpath, subid{ii}, subdir, maskname);
    % Load mask with MRIread
    headeronly = 0;
    permuteflag = 1;
    vol = MRIread(maskpath,headeronly, permuteflag);
    mask = vol.vol;
    %% Save the mask matrix as a .MAT
    % Names of output files
    maskout = fullfile(dpath, subid{ii}, subdir, 'maskec.mat');
    % Save as .MAT
    save(maskout, 'mask', '-v7.3');
end