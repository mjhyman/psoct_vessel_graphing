%% LSM Vasculature Analyses
%{
- Use the vessel graphing
- Analyze geometry
%}


%% Import LSM tif
clear; clc; close all
% Data paths
dpath = '/projectnb/npbssmic/ns/LSM/';
fname = '2068_Ilastik_Simple_Segmentation_OME_Resized.tif';

% Add upper-level path
top_path = '/projectnb/npbssmic/s/mhyman/bunpc_psoct_vessel_graphing/psoct_vessel_graphing/';
addpath(genpath(top_path));

% Import TIF image as binary
fprintf('\nImporting TIF file\n')
lsm = imbinarize(tiffreadVolume(fullfile(dpath, fname)));

%% Image Quality Assurance
% Remove small, disjoint sections
% Define a threshold for removing small disjoint sections
nvox = 50;
% Remove small disjoint sections from the binary image
lsm_mod = bwareaopen(lsm, nvox);

%% Separate into large, medium, small networks

% Find connected components
cc = bwconncomp(lsm_mod);
numPixels = cellfun(@numel, cc.PixelIdxList);

% 1. Find the largest connected component
fprintf('\nDefining largest network\n')
[~, maxIdx] = max(numPixels);
net_large = false(size(lsm_mod));
net_large(cc.PixelIdxList{maxIdx}) = true;

% 2. Define Thresholds using prctile (Standard MATLAB)
% Adjust these percentiles based on your specific distribution
nlow = prctile(numPixels, 30); 
nhigh = prctile(numPixels, 60);

% 3. Extract Medium Networks
% Find indices where the component size is within the range
fprintf('\nDefining medium network\n')
med_indices = find(numPixels > nlow & numPixels <= nhigh);
net_med = false(size(lsm_mod));
if ~isempty(med_indices)
    % Combine all pixel indices for these components into one vector
    net_med(vertcat(cc.PixelIdxList{med_indices})) = true;
end

% 4. Extract Small Networks
% Find indices where component size is below the low threshold
% small_indices = find(numPixels <= nlow);
% net_small = false(size(lsm_mod));
% if ~isempty(small_indices)
%     % Combine all pixel indices for these components
%     net_small(vertcat(cc.PixelIdxList{small_indices})) = true;
% end

%% Skeletonize + convert to graph
% Minimum branch length = 10 voxels
min_branch_len = 10;

% Initialize struct for storing nodes/edges/connections
graphs = struct();

% Largest network
fprintf('\nGraphing Largest Network\n')
[nodes,edges,conns] = vol_to_graph(net_large,min_branch_len);
graphs.large.nodes = nodes;
graphs.large.edges = edges;
graphs.large.conns = conns;

% Medium networks
fprintf('\nGraphing Medium Networks\n')
[nodes,edges,conns] = vol_to_graph(net_med,min_branch_len);
graphs.med.nodes = nodes;
graphs.med.edges = edges;
graphs.med.conns = conns;

%% Analyze the graphs

%% Export largest network into 3D volume
% Downsample as necessary for visualization