function [ labeled_elements, gridbmap, idx_neighbors, angles ] = obtain_contour_elements( seg, vectors )
% [ labeled_elements, gridbmap, idx_neighbors, angles ] = obtain_contour_elements( seg, vectors )
% ------------------------------------------------------------------------
% Obtains the elements that form the contours of the segmentation.
%
% INPUT
%	seg                    Input segmentation.
%	vectors                Contour vectors resulting from the function 
%                          create_finest_partition (performed during 
%                          the ucm creation).
%
% OUTPUT
%	labeled_elements	   Indexed contour elements.
%   gridbmap               Position of the contour elements.
%   idx_neighbors          Regions that share each contour element.
%   angles                 Direction of each contour element.
%
%
% ------------------------------------------------------------------------
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
%
%  David Varas <david.varas@upc.edu>
%  May 2014
% ------------------------------------------------------------------------

[gridbmap, idx_neighbors, junctions] = seg2gridbmap(seg,0);

[gridbmap, idx_neighbors_matrix_min, idx_neighbors_matrix_max, angles ] = select_gridbmap(gridbmap, idx_neighbors.matrix_min, idx_neighbors.matrix_max, vectors);

gridbmap = gridbmap > 0;

idx_neighbors.matrix_min = idx_neighbors_matrix_min;
idx_neighbors.matrix_max = idx_neighbors_matrix_max;

labeled_elements = label_contour_elements(double(gridbmap));


end

