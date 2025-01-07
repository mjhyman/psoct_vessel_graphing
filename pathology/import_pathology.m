function [stain, mask, res] = import_pathology(dpath, stain_name, mask_name)
%IMPORT_PATHOLOGY read the TIF file for the pathology (A-beta or p-tau)
% Retrieve TIFF metadata for the staining & compute the resolution of each
% pixel. Then, load the TIFF and the respective tissue mask. Rescale the
% staining to [0,1] to standardize across all images. Then, apply the
% tissue mask.
% INPUTS:
%   dpath (string): directory path to the staining TIFF file
%   stain_name (string): filename of staining TIFF file
%   mask_name (string): filename of tissue mask TIFF file
% OUTPUTS:
%   stain (double matrix): masked staining, scaled between [0,1]
%   mask (logical matrix): tissue mask for the pathology slide
%   bg (double matrix): background patch for normalizing the pathology.
%                       This is a region of the stain from the white matter
%                       that does not include any plaques.
%   res (double array): resolution of a pixel [x,y]

%%% Disable Tiff read warnings for unknown tags
id = 'imageio:tiffmexutils:libtiffWarning';
warning('off',id);

%%% Retrieve metadata and calculate resolution
% Read staining TIFF info
stain_file = append(stain_name, '.tif');
stain_file = fullfile(dpath, stain_file);
stain_metadata = imfinfo(stain_file);
% X and Y resolution (pixels / micron)
xres = stain_metadata.XResolution;
yres = stain_metadata.YResolution;
% Invert resolution (microns / pixel)
xres = 1./xres;
yres = 1./yres;
res = [xres, yres];

%%% Load TIFF file
% Read staining TIFF (RGB) & convert to grayscale
stain = Tiff(stain_file,'r');
stain = read(stain);
stain = rescale(stain, 0, 1);
% stain = rgb2gray(stain);
% Load tissue mask
mask_fname = append(mask_name, '.tif');
mask_file = fullfile(dpath, mask_fname);
mask = Tiff(mask_file,'r');
mask = logical(read(mask));
end
