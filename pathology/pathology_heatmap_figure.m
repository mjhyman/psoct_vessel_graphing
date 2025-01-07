%% Generate Heatmap of Cropped p-tau
% Overview:
%   - generate heatmap of a cropped image of the p-tau stain from subject
%   AD_21354. This is for a figure in the AD/CTE paper.

%% Add top-level directory of code repository to path
clear; clc; close all;
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Add top-level directory to path to find functions
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));
% Set maximum number of threads equal to number of threads for script
NSLOTS = str2num(getenv('NSLOTS'));
maxNumCompThreads(NSLOTS);

%% Initialize directories, filenames, parameters
% Path to pathology image
pathology_path = ['/projectnb/npbssmic/pantinew/tau_amyloid_images/AT8' ...
                  '/AD 21354 AT8/'];
fname = 'crop_p-tau_threshold.tif';

%% Call function to create heatmap
% Import pathology cropped
stain = imcomplement(TIFF2MAT(fullfile(pathology_path, fname)));

% Rescale stain
stain = rescale(stain,0,1);
% Create mask to use entirety of image
mask = ones(size(stain));

% Declare pixel size of pathology image (microns)
pix = [0.2194598,0.2194598];
% Declare size of averaging box in heatmap (microns)
square_side = 10;

%% Function to create heatmap

%%% Create heatmap
% Create heatmap
hm = pathology_heatmap(pix, square_side, stain, mask);
% Plot heatmap
figure; imagesc(hm);

%%% Adjust colors
% Initialize the colormap limits
cmap = jet(256);
clim(gca, [0, 1]);
colormap(cmap);
% Add colorbar
colorbar;

%%% Add scale bar
% Set length (microns)
scale_bar = 250;
% Calculate number of pixels for scale bar
npix = scale_bar / pix(1);
% Add scale bar
yoff = size(hm,1) - 200;
xoff = size(hm,1) - npix - 25;
hold on;
plot([xoff; xoff + npix],[yoff; yoff],'-w','LineWidth',10);
% Set font size of scale bar
set(gca,'FontSize',24)
% Remove tick labels
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);


%%% Save the output
fout = 'crop_p-tau_threshold_heatmap_scale_250um';
fout = fullfile(pathology_path,fout);
saveas(gca,fout,'png')












