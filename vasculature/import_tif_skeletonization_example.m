%% Demo script for 3D skeletonization
% This script performs the following:
% - import a .TIF of a 3D volume
% - convert the volume into a logical
% - run the Matlab built-in function "bwskel" to identify the centerline of
%       each group of voxels.
% - Visualize the output

%% Import volume & convert to logical (binary)

% Import volume
filepath = 'path_to_file';
filename = 'filename.tif';
dtype = 'uint8';
vol = logical(import_tif(filepath,filename));

%% Create skeleton centerline w/ bwskel & visualize

% Function to skeletonize the logical volume
skel = bwskel(vol);
% Figure title for visualization
fig_title = 'overlay skeleton (red) and segmentation (gray)';
% Overlay of skeleton and segmentation
[~,~] = overlay_seg_skel(seg,skel,fig_title);

%% Function to import 3D .TIF into matlab
function image_stack = import_tif(filepath, filename,dtype)

% Create full filepath to 3D volume
full_file = fullfile(filepath, filename);
% Create TIFF object
t = Tiff(full_file, 'r');
info = imfinfo(filename);
num_slices = numel(info);
% Initialize imaging stack
image_stack = zeros(info(1).Height, info(1).Width, num_slices, dtype);
% Iterate over depths
for k = 1:num_slices
    t.setDirectory(k);
    image_stack(:, :, k) = t.read();
end
t.close();
end

%% Function to overlay skeleton and segmentation
function [seg, skel] = overlay_seg_skel(seg, skel, fig_title)
% overlay the skeleton and segmentation
%   This code overlays the 3D figure of the segmentation + skeleton
% INPUTS:
%   seg (matrix): segmentation (binary)
%   skel (matrix): skeleton of graph (binary)
%   fig_title (string): name of figure
% OUTPUTS:
%   displays a 3D volume of segmentation and skeleton

%% Verify dimensions of skeleton and segmentation matrices are equivalent
% Skeleton dimensions
skel_sz = size(skel);
% Segmentation dimensions
seg_sz = size(seg);
% Assert the two dimensions match, otherwise raise an error
assert(all(skel_sz == seg_sz), ['The dimensions of the skeleton and ' ...
    'the segmentation do not match.']);

%% If skeleton matrix dim exceeds 2048 -> downsample skeleton & segmentation
% tmp will equal 1 if any of the dimensions exceed limits of "viewer3D"
% function (which is 2048)
maxdim = 2048;
sz = skel_sz;
tmp = sz > maxdim;
if any(tmp)
    %%% Calculate scaling factor
    % Find largest dimension exceeding 2048
    m = max(sz(tmp));
    % Calculate smallest scaling factor to reach 2048
    scalar = ceil(m ./ maxdim);
    
    %%% Perform scaling according to scaling factor.
    % Maintain z-axis dimension, so long as it doesn't exceed 2048
    if sz(3) < 2048
        %%% Segmentation reshape and squeeze (fresh lemonade)
        seg_re = reshape(seg, scalar, sz(1)/scalar, scalar,...
                            sz(2)/scalar, 1, sz(3));
        seg = squeeze(sum(seg_re, [1,3,5])) / scalar.^2;
        %%% Skeleton reshape and squeeze (fresh lemonade)
        skel_re = reshape(skel, scalar, sz(1)/scalar, scalar,...
                            sz(2)/scalar, 1, sz(3));
        skel = squeeze(sum(skel_re, [1,3,5])) / scalar.^2;
    % Downsample z-axis dimension, if its => 2048
    else
        %%% Skeleton reshape and squeeze (fresh lemonade)
        seg_re = reshape(seg, scalar, sz(1)/scalar, scalar,...
                            sz(2)/scalar, scalar, sz(3)/scalar);        
        seg = squeeze(sum(seg_re, [1,3,5])) / scalar.^3;
        %%% Skeleton reshape and squeeze (fresh lemonade)
        skel_re = reshape(skel, scalar, sz(1)/scalar, scalar,...
                            sz(2)/scalar, scalar, sz(3)/scalar);        
        skel = squeeze(sum(skel_re, [1,3,5])) / scalar.^3;
    end
end

%% Initialize the 3D figure and display

%%% Initialization parameters
view_panel = uifigure(figure,'Name',fig_title); close;
v = viewer3d(view_panel);
v.BackgroundColor = 'w';
v.BackgroundGradient = 'off';

%%% Display volume of skeleton (from graph)
h = volshow(skel,'Parent',v);
h.Parent.BackgroundColor = 'w';
% Make skeleton red
skelcolor = repmat([1, 0, 0], [256,1]);
h.Colormap = skelcolor;

%%% Overlay the segmentation
h.OverlayData = seg;
h.OverlayAlphamap = 0.1;
% Make overlay gray
gray = repmat([0.5, 0.5, 0.5], [256,1]);
h.OverlayColormap = gray;

end