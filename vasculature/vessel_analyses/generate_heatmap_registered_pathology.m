%% Heatmap of the registered pathoogy image
% Overview:
%   - import the registered pathology heatmaps
%   - determine the dimensions of each heatmap
%   - use plotting function that scales each pathology heatmap accordingly

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
% Remove the two sub folders to reach parent
% (psoct_human_brain\vasculature\vesSegment)
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));
% Set maximum number of threads equal to number of threads for script
NSLOTS = str2num(getenv('NSLOTS'));
maxNumCompThreads(NSLOTS);

%% Initialize directories, filenames, parameters
%%% Directories 
% Metrics output path
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];
% Registered staining heatmaps directory
path_reg = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/heatmaps/' ...
    'heatmaps_pathology_registration/'];

%%% Vasculature heatmap matrix filename
% ROI cube side (microns)
cube_side = 1000;
% Load according to size of ROI cube
hm_fname = append('heatmap_ab_ptau_',num2str(cube_side),'.mat');
% Load the vascular heatmap matrix
hm = load(fullfile(mpath,hm_fname));
hm = hm.heatmap;
subid = fields(hm);

%%% Initialize filenames of registered pathology / heatmaps.
% Some of the subjects had the patholoy registered with the automated
% "demons" algorithm and other were registered with the Maritos manual
% landmark software. This section will specify this delineation.
%
% Subjects with heatmaps in TIF
ftif_ab_path = 'Ab_registered.tif';
ftif_ab_mask = 'Ab_mask_registered.tif';
ftif_pt_path = 'AT8_registered.tif';
ftif_pt_mask = 'AT8_mask_registered.tif';
subs_tif = {'AD_20832','AD_20969','CTE_6489','CTE_6912',...
    'NC_6839','NC_21499'};
% Subjects with heatmaps in .MAT (the actual filenames have the subject ID
% appended as a prefix (i.e. [sub_id]_[file_name])
fmat_ab_path = 'Ab_path_registered.mat';
fmat_ab_mask = 'Ab_path_mask_registered.mat';
fmat_pt_path = 'AT8_path_registered.mat';
fmat_pt_mask = 'AT8_path_mask_registered.mat';
subs_mat = {'AD_10382','AD_21354','AD_21424','CTE_7019','CTE_7126',...
    'NC_8095'};

%%% Subvolume parameters for analyses
% Isotropic cube length (microns)
cube_side = 1000;
% Size of each OCT voxel (microns)
vox = [12, 12, 15];
% Whether to plot non-normalized heatmaps for each depth
viz_individual = false;
% Pathology isotropic cube side length (microns) 
path_cube_side = 204;
% Pixel size of pathology
res = [10.9731, 10.9731];

%%% Compute number of voxels in x,y,z for each cube
n_x = floor(cube_side ./ vox(1));
n_y = floor(cube_side ./ vox(2));
n_z = floor(cube_side ./ vox(3));
% Calculate the size of each cube (in voxels)
cube_vol_vox = n_x * n_y * n_z;
% Calculate the size of each cube (in cubic microns)
cube_vol_um = cube_vol_vox * vox(1) * vox(2) * vox(3);

%%% Index for structs
% Pathology
pidx = {'ab','pt'};
% masks (gray matter and white matter)
midx = {'gm','wm'};

%%% Struct for storing pairs of values (heatmap, pathology)
pairs = struct();

%%% label cell array
xlabels = {'Volume Fraction (unitless)','Branch Density (mm^-^3)',...
    'Length Density (mm / mm^3)','Tortuosity (unitless)'};
ylabels = {'[A-beta]','[p-tau]'};

%% Calculate maximum heatmap dimensions
% This will be used for scaling each figure so that all heatmaps are on
% the same scale. This simplifies figure generation.

% initialize the maximum values for the x,y dimensions
ymax = 1; xmax = 1;

% Iterate over each subject
for ii = 1:length(subid)
    %%% Load the heatmap for the subject
    sub = subid{ii};
    % Identify the maximum dimensions of x,y
    ymax = max([ymax,size(hm.(sub).mask,1)]);
    xmax = max([xmax,size(hm.(sub).mask,2)]);
end

%% Figure of pathology heatmap
% Load the vascular heatmap and pathology heatmaps, recreate the ROIs from
% the vascular heatmaps, generate an x-y pair for each ROI (x = vascular
% ROI, y = pathology ROI), perform Spearman's on the matrix pairs for each
% subject and stain.

% Iterate over each subject
for ii = 1:length(subid)
    %% Load the registered pathology and retrieve vascular heatmaps
    % Retrieve the current subject ID
    sub = subid{ii};
    % Determine which registration method was used on this subject 
    if any(ismember(subs_tif,sub))
        % TIF: Load the A-beta and p-tau pathology heatmaps
        ab = TIFF2MAT(fullfile(path_reg,sub,ftif_ab_path));
        ab_mask = TIFF2MAT(fullfile(path_reg,sub,ftif_ab_mask));
        pt = TIFF2MAT(fullfile(path_reg,sub,ftif_pt_path));
        pt_mask = TIFF2MAT(fullfile(path_reg,sub,ftif_pt_mask));
        % Convert from RGB to grayscale (all depths are equivalent
        ab = ab(:,:,1); ab_mask = logical(ab_mask(:,:,1));
        pt = pt(:,:,1); pt_mask = logical(pt_mask(:,:,1));
    else
        % MAT: Load the A-beta and p-tau pathology heatmaps
        ab = load(fullfile(path_reg,sub,append(sub,'_',fmat_ab_path)));
        ab_mask = load(fullfile(path_reg,sub,append(sub,'_',fmat_ab_mask)));
        pt = load(fullfile(path_reg,sub,append(sub,'_',fmat_pt_path)));
        pt_mask = load(fullfile(path_reg,sub,append(sub,'_',fmat_pt_mask)));
        % Load the fields of each struct
        ab = ab.path_registered; pt = pt.path_registered;
        ab_mask = logical(ab_mask.path_mask_registered);
        pt_mask = logical(pt_mask.path_mask_registered);
    end

    %% Add to struct
    % Retrieve the vascular heatmaps/masks from the "hm" struct
    vasc = hm.(sub);
    masks = struct();
    masks.gm = logical(vasc.mask_gm);
    masks.wm = logical(vasc.mask_wm);
    hm_vf = vasc.vf;
    hm_bd = vasc.bd;
    hm_ld = vasc.ld;
    hm_tr = vasc.tort;
    hm_dm = vasc.diam;

    % Combine the masks and pathology into struct
    path = struct();
    path.ab.mask = ab_mask;
    path.pt.mask = pt_mask;
    path.ab.hm = ab;
    path.pt.hm = pt;

    %% Pathology heatmap
    %%% Compute heatmap from stain
    % Rescale a-beta and p-tau
    ab_re = rescale(ab,'InputMin',0,'InputMax',2^8);
    pt_re = rescale(pt,'InputMin',0,'InputMax',2^8);
    % Compute heatmaps for A-beta and p-tau
    [hm_ab] = pathology_heatmap(res, path_cube_side, ab_re, ab_mask);
    [hm_pt] = pathology_heatmap(res, path_cube_side, pt_re, pt_mask);
    
    %%% Rotate / invert subjects to align gyri at top of figure
    % These rotations match those in the vascular heatmaps
    if strcmp(sub,'AD_20969')
        hm_ab = flip(permute(hm_ab,[2,1,3]),2);
        ab_mask = flip(permute(ab_mask,[2,1,3]),2);
        hm_pt = flip(permute(hm_pt,[2,1,3]),2);
        pt_mask = flip(permute(pt_mask,[2,1,3]),2);
    end
    if strcmp(sub,'CTE_6912')||strcmp(sub,'NC_6974')||strcmp(sub,'NC_21499')
        hm_ab = flip(hm_ab,1);
        ab_mask = flip(ab_mask,1);
        hm_pt = flip(hm_pt,1);
        pt_mask = flip(pt_mask,1);
    end
    if strcmp(sub,'CTE_6489')
        hm_ab = imrotate(hm_ab,-15);
        ab_mask = imrotate(ab_mask,-15);
        hm_pt = imrotate(hm_pt,-15);
        pt_mask = imrotate(pt_mask,-15);
    end

    %%% Plot/save heatmap
    % Output path
    fpath = fullfile(path_reg,sub);
    % Figure Names
    ab_fname = append(sub,'_Ab_registered_stain_heatmap_',...
        num2str(path_cube_side),'.png');
    pt_fname = append(sub,'_AT8_registered_stain_heatmap_',...
        num2str(path_cube_side),'.png');
    ab_title = append(sub,' A-beta');
    pt_title = append(sub,' p-tau');
    % Plot heatmaps
    plot_save_heatmap(1,hm_ab,0,[0,1],[ymax,xmax],ab_mask,ab_title,...
        '% Area',fpath,ab_fname)
    plot_save_heatmap(1,hm_pt,0,[0,1],[ymax,xmax],pt_mask,pt_title,...
        '% Area',fpath,pt_fname)

end



%% Plot and save the heat maps
function plot_save_heatmap(Ndepths, heatmaps, flip_cbar, colorbar_range,...
    max_dim, masks, tstr, cbar_label, dpath, fname)
% PLOT_SAVE_HEATMAP use imagesc and set background = 0
% INPUT
%   Ndepths (int): number of depths in z dimension
%   heatmaps (double matrix): heatmaps of vascular metric
%   flip_cbar (logical): reverse the direction of the colorbar
%   colorbar_range (double array): [min, max]
%   max_dim (double array): [ymax, xmax]
%   masks (double): tissue mask (1=tissue, 0=other)
%   tstr (string): figure title
%   cbar_label (string): colorbar label
%   dpath (string): data directory path
%   fname (string): name of figure to save

%%% Set the number of depths to iterate for each heatmap
% If the number of depths is not specified, then set it equal to the number
% of z dimensions.
if isempty(Ndepths)
    Ndepths = size(heatmaps,3);
end
% Set fontsize for the heatmap figure
fontsize = 20;
% Set the maximum dimensions for x-y axes
ymax = max_dim(1);
xmax = max_dim(2);

%%% Iterate over frames in z dimension
for d = 1:Ndepths
    %%% Heatmap of j_th frame from the length density
    % If there are multiple heatmaps in the matrix
    if size(heatmaps,3) > 1
        heatmap = heatmaps(:,:,d);
    % Here it is just a single frame of a heatmap
    else
        heatmap = heatmaps;
    end
    % Initialize heatmap
    h = imagesc(heatmap);

    %%% Scale the x- and y-axes according to the max dimensions
    % Calculate the ratio of the current axes to the max axes
    yratio = size(heatmap,1) ./ ymax;
    xratio = size(heatmap,2) ./ xmax;
    set(gcf,'Units','Normalized','OuterPosition',[0,0,xratio,yratio])

    %%% Initialize colorbar
    % If the colorbar_range is passed in, then extract min & max
    if ~isempty(colorbar_range)
        cmap_min = colorbar_range(1);
        cmap_max = colorbar_range(2);
    % Otherwise, set limits from the current heatmap
    else
        % Min = Lowest value also greater than zero
        cmap_min = min(heatmap(heatmap(:)>0));
        % Max = Find the 95th percentile for upper limit
        cmap_max = prctile(heatmap(:), 95);
    end
    % Initialize the colormap limits
    cmap = jet(256);
    clim(gca, [cmap_min, cmap_max]);
    % Initialize colormap and colorbar
    if flip_cbar
        colormap(flipud(cmap));
    else
        colormap(cmap);
    end
    c = colorbar;

    %%% Apply tissue mask to the heatmap to remove background
    alpha_mask = double(masks(:,:,d));
    set(h, 'AlphaData', alpha_mask);

    %%% Configure figure parameters
    % Update title string with specific pathology
    pathology = {'A-Beta','p-tau'};
    if size(heatmaps,3) > 1
        title_str = append(tstr, ' (', pathology{d},' depth)');
    else
        title_str = tstr;
    end
    title(title_str,'Interpreter','none');
    set(gca, 'FontSize', fontsize);
    % Label the colorbar    
    c.Label.String = cbar_label;
    % Offset colorbar label to the right of colorbar
    c.Label.Position = [10 (cmap_max - (cmap_max-cmap_min)/2)];
    c.Label.Rotation = 270;
    % Increase fontsize of colorbar
    c.FontSize = 20;
    % Remove x and y tick labels
    set(gca,'Yticklabel',[]);
    set(gca,'Xticklabel',[]);
    
    %%% Save figure as PNG
    % If there are multiple heatmaps in the matrix, save vascular heatmap
    if size(heatmaps,3) > 1
        fout = append(fname, '_', pathology{d});
    % Otherwise, save the pathology heatmap
    else
        fout = fname;
    end
    % If the colorbar is reversed then add suffix to filename
    if flip_cbar
        fout = append(fout, '_flip_cbar');
    end
    % Save figure as PNG
    fout = fullfile(dpath, fout);
    pause(1)
    box off;
    saveas(gca, fout,'png');
    close;
end
end