function [ image_i, image_j ] = create_image_segmentations( seg_i, seg_j )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

regions_i = max(max(seg_i));
regions_j = max(max(seg_j));

max_regions = max(regions_i,regions_j);

tmp = jet(max_regions);
color = uint8(255*tmp(randperm(max_regions,max_regions),:));

[s_x, s_y] = size(seg_i);
for i=1:s_x
    for j=1:s_y
        image_i(i,j,1) = color(seg_i(i,j),1);
        image_i(i,j,2) = color(seg_i(i,j),2);
        image_i(i,j,3) = color(seg_i(i,j),3);
        
        
        image_j(i,j,1) = color(seg_j(i,j),1);
        image_j(i,j,2) = color(seg_j(i,j),2);
        image_j(i,j,3) = color(seg_j(i,j),3);
    end
end
end

