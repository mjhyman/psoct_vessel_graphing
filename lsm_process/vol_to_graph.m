function [nodes, edges, adj_list] = vol_to_graph(vol, min_branch)
    % vol_to_graph Convert 3D binary into a node-edge graph
    % OUTPUTS:
    %   [nodes]: N x 3 matrix of [y, x, z] coordinates
    %   [edges]: Cell array containing the [y, x, z] voxels for each segment
    %   [adj_list]: M x 2 matrix of node index pairs (connections)
    % INPUTS:
    %   vol: binary volume
    %   min_branch (uint): minimum branch length in voxels

    %% 0. Skeletonize the binary volume
    skel = bwskel(vol,'MinBranchLength',min_branch); 

    %% 1. Identify Raw Node Voxels
    kernel = ones(3,3,3);
    kernel(2,2,2) = 0;
    % Count neighbors for every skeleton pixel
    neighbor_count = convn(double(skel), kernel, 'same') .* double(skel);
    
    % Junctions (>2 neighbors) and Endpoints (1 neighbor)
    raw_node_mask = (neighbor_count > 2) | (neighbor_count == 1);
    
    %% 2. Pro-Tip: Cluster Junction Voxels
    % bwskel often creates clusters of junction pixels at one intersection.
    % We group these and take the centroid to create ONE node per intersection.
    cc_nodes = bwconncomp(raw_node_mask);
    nodes = zeros(cc_nodes.NumObjects, 3);
    
    % Map every skeleton node-pixel to its new unique "Cluster ID"
    pixel_to_node_map = zeros(size(skel));
    
    for i = 1:cc_nodes.NumObjects
        % Get coordinates of pixels in this cluster
        [y, x, z] = ind2sub(size(skel), cc_nodes.PixelIdxList{i});
        % Store the centroid as the official node location
        nodes(i, :) = [mean(y), mean(x), mean(z)];
        % Mark these pixels in the map
        pixel_to_node_map(cc_nodes.PixelIdxList{i}) = i;
    end

    %% 3. Extract Edges (Segments)
    % Remove all raw node pixels to isolate the segments
    edge_skeleton = skel;
    edge_skeleton(raw_node_mask) = 0;
    
    cc_edges = bwconncomp(edge_skeleton);
    num_edges = cc_edges.NumObjects;
    edges = cell(num_edges, 1);
    adj_list = [];

    %% 4. Trace Connectivity
    % Use a 3x3x3 dilation to see which nodes each edge segment touches
    se = ones(3,3,3);
    
    for i = 1:num_edges
        % Store edge voxel coordinates
        [ey, ex, ez] = ind2sub(size(skel), cc_edges.PixelIdxList{i});
        edges{i} = [ey, ex, ez];
        
        % Create a temporary mask for this specific edge
        temp_edge_mask = false(size(skel));
        temp_edge_mask(cc_edges.PixelIdxList{i}) = true;
        
        % Dilate the edge to find adjacent node clusters
        dilated_edge = imdilate(temp_edge_mask, se);
        intersecting_pixels = find(raw_node_mask & dilated_edge);
        
        % Identify unique Node IDs (from our map) this edge touches
        connected_node_ids = unique(pixel_to_node_map(intersecting_pixels));
        % Clean up (remove zeros if any)
        connected_node_ids(connected_node_ids == 0) = [];
        
        % If the edge connects nodes, add to adjacency list
        if length(connected_node_ids) >= 2
            % Usually 2 nodes, but can be more in rare high-complexity clusters
            % We add the pair [StartNode, EndNode]
            adj_list = [adj_list; connected_node_ids(1), connected_node_ids(2)];
        elseif length(connected_node_ids) == 1
            % This is a "dangling" edge or a loop back to the same node
            adj_list = [adj_list; connected_node_ids(1), connected_node_ids(1)];
        end
    end
end