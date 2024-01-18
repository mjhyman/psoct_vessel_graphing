function [volm, t] = create_mask_v2(vol, debug, NSLOTS)
%CREATE_MASK create/apply mask for each slice
% Create boundary mask, imerode (diamond), output the mask
% INPUTS:
%   vol (double matrix): psoct c-scan of tissue volume
% OUTPUTS:
%   volm (uint16): stack of mask for each slice
%   t (double): array of times (seconds) to run each active contour
% Mathworks References:
%   Multithresholding: Otsu, N., "A Threshold Selection Method from Gray-
%                   Level Histograms," IEEE Transactions on Systems, Man,
%                   and Cybernetics, Vol. 9, No. 1, 1979, pp. 62-66.
%   Active Contour: T. F. Chan, L. A. Vese, Active contours without edges.
%                   IEEE Transactions on Image Processing, Volume 10, Issue
%                   2, pp. 266-277, 2001.

%% Initialize Parameters and Variables
% Initialize matrix for storing masked volume
volm = uint16(zeros(size(vol)));
% Create structuring element for eroding the mask
se = strel('disk',3);
% Number of iterations for active contour edge finding
nactive = 500;
% Define array for storing timer data to track perforamnce of activecont.
t = zeros(1, size(vol, 3));
% Set the range of min/max object size (connected pixels) to retain.
min_size = 100;
max_size = size(vol,1) * size(vol,2);
range = [min_size, max_size];

%% Iterate through volume stack and find boundary mask for each slice
for ii = 1:size(vol, 3)
    %% Initialize Slice from volume + Debugging info
    s = vol(:,:,ii);
    %%% Debugging information
    fspec = 'Starting slice %i\n';
    fprintf(fspec, ii)
    
    %% Counter for b-scan within physical stack
    % Variable for number of images per physical stack
    phy = 11;
    % Array for index within physical stack
    stackidx = mod(ii,phy+1);
    % If the stack index is 0 or between 4-10, then do not compute a mask.
    % This corresponds to an image index of 4-11.
    if (1 <= stackidx) && (stackidx <= 3)
        %%% Initialize mask -- multithreshold slice to find agarose pixels
        lvl = multithresh(s,5);
        % Discard the darkest portions of image (lowest threshold)
        lvl = lvl(2:end);
        % Quantize image with thresholds
        maskth = imquantize(s, lvl);
        % Set (quanta > 1) = 1 and (quanta == 1) = 0
        maskth(maskth==1) = 0;
        maskth(maskth>1) = 1;
        % Fill mask
        maskfill = imfill(maskth);
        % Remove islands
        mask0 = bwareafilt(logical(maskfill), range,4);
        %%% Set the method for active contour
        % The "edge" method tends to contract to meet the edge
        ac_meth = 'edge';
    else
        %%% Initialize mask by registering mask from third slice
        % Select the third slice in the physical slice
        idx = ii - (stackidx - 3);
        %{
        moving = vol(:,:,idx);
        % Register slices (fixed = current slice, moving = third slice)
        tform = imregcorr(moving, s);
        % Create image dimension reference from slice
        rfixed = imref2d(size(moving));
        % Apply registration to mask from third slice in sequence
        mask0 = imwarp(volm(:,:,idx), tform, "Outputview", rfixed);
        %}
        % Use the third mask in physical slice series
        mask0 = volm(:,:,idx);

        %%% Set the method for active contour
        % The Chan-Vese method can expand or contract the mask to meet the
        % edges, whereas the "edge" method can only contract.
        ac_meth = 'edge';
    end

    %% Iteratively find edge, erode border, remove islands
    %%% Active contour - find tissue border
    tic;
    ac = activecontour(s, mask0, nactive, ac_meth);
    t(ii) = toc;
    % Output the time for each iteration
    fspec = 'Iteration %i took %i seconds to run active contour\n';
    fprintf(fspec, ii, t(ii))
    
    %%% Erode the border + Remove islands of pixels (disjoint pixels)
    bw = imerode(ac, se);
    volm(:,:,ii) = bwareafilt(bw, range);

    %% Debugging figures
    if debug
        % If image 1-3 in stack, then plot these
        tmp = exist('maskth','var');
        if tmp
            figure; imshow(maskth, []); title('Multithresholded')
            figure; imshow(maskfill, []); title('Multithresholded + Filled')
        end
        % Clear these variables before next iteration
        clear maskth
        % Initial mask to feed into active contour
        figure; imshow(mask0, []);
        title('Initial Mask - ',['Slice ', num2str(ii)])
        % Active contour
        figure; imshow(ac, []);
        title({'Active contour',['Slice ', num2str(ii)]})
        % Active contour + borders eroded
        figure; imshow(bw, []);
        title({'Active contour, borders eroded, islands removed',...
            ['Slice ', num2str(ii)]})
        % Active contour + borders eroded + islands removed
        figure; imshow(volm(:,:,ii), []);
        title({'Active contour, border eroded, islands removed',...
            ['Slice ', num2str(ii)]})
    end

end

end
