function [hm] = pathology_heatmap(pix, square_side, stain, mask)
%PATHOLOGY_HEATMAP generate heatmap with pathology
% Divide the pathology stain into isotropic squares, calculate the density
% of each square, add the density to a matrix.
%
% INPUTS:
%   pix (double array): pixel size (x,y) (microns)
%   square_side (uint): heatmap square dimension (microns)
%   stain (uint8 matrix): pathology stain
%   th (double): minimum threshold for segmenting AB or AT8 (after taking
%               complement of the staining)
%   mask (logical matrix): tissue mask for the pathology stain
% OUTPUTS:
%   hm (double matrix): masked heatmap of the pathology stain

%%% Verify that the stained image is scaled between [0,1]
% The inversion assumes that the image is scaled between [0,1], so this
% will return the incorrect value otherwise.
r = range(stain);
assert(0<=min(r), 'The stained image is not scaled between [0,1]');
assert(max(r)<=1, 'The stained image is not scaled between [0,1]');

%%% Translate cube size to pixel size
% Number of pixels in each dimension to create square
nx = floor(square_side ./ pix(1));
ny = floor(square_side ./ pix(2));

%%% Generate heat map
% Initialize heatmap matrices for A-beta and p-tau
hm = zeros(size(stain,1), size(stain,2));
% Iterate over rows
for x = 1:nx:size(stain,1)
    % Iterate over columns
    for y = 1:ny:size(stain,2)
        %%% Crop segmentation into isotropic square
        % Initialize end indices for each axis
        xf = x + nx - 1;
        yf = y + ny - 1;
        % Take minimum of matrix dimensions and end indices
        xf = min(xf, size(stain,1));
        yf = min(yf, size(stain,2));
        % Take cube from segmentation
        path_square = stain((x:xf), (y:yf));
        % Take cube from mask
        mask_square = mask((x:xf), (y:yf));
        % Remove plaque pixels outside of mask boundary
        path_square(~mask_square) = 0;

        %%% Calculate propotion of plaque pixels in subvolume
        % Calculate the total nonzero area of the mask subvolume
        m_area = sum(mask_square(:));
        % Calculate the total nonzero area of the plaque subvolume
        p_area = sum(path_square(:));
        % Calculate the ratio of plaque : mask pixels
        hm((x:xf), (y:yf)) = p_area ./ m_area;
    end
end

% Mask the pathology heatmap with tissue mask
hm(~mask) = 0;
end