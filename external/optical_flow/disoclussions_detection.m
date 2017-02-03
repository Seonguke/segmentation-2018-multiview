function [ disocluded_regions ] = disoclussions_detection( I_i, I_j, seg_i, OF_regions_x, OF_regions_y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

regions_i = max(max(seg_));
[s_x, s_y] = size(seg_i);

disocluded_regions = zeros(1, regions_i);

for i=1:regions_i
    region_mask_i = (seg_i==i);
    region_mask_j = zeros(s_x, s_y);
    
    delta_x = round(OF_regions_x(i));
    delta_y = round(OF_regions_y(i));
    
    for u=1:s_x
        for v=1:s_y
            if(region_mask_j(u,v)) 
                if((u-delta_y)>0 && (u+delta_y)<=s_y && (v+delta_x)>0 && (v+delta_x)<=s_x)
                    region_mask_j(u+delta_y,v+delta_x) = 1;
                end
            end
        end
    end
    
    %% Extract red pixels
    image_mask = zeros(s_x, s_y, 3);
    image_mask(:,:,1) = region_mask_i;
    red_pixels_i = I_i(image_mask==1);
    
    image_mask(:,:,1) = region_mask_j;
    red_pixels_j = I_j(image_mask==1);
    
    %% Extract green pixels
    image_mask = zeros(s_x, s_y, 3);
    image_mask(:,:,2) = region_mask_i;
    green_pixels_i = I_i(image_mask==1);
    
    image_mask(:,:,2) = region_mask_j;
    green_pixels_j = I_j(image_mask==1);

    %% Extract blue pixels
    image_mask = zeros(s_x, s_y, 3);
    image_mask(:,:,3) = region_mask_i;
    blue_pixels_i = I_i(image_mask==1);
    
    image_mask(:,:,3) = region_mask_j;
    blue_pixels_j = I_j(image_mask==1);
    
    %% Compute means
    diff_red   = abs(mean(red_pixels_i) - mean(red_pixels_j));
    diff_green = abs(mean(green_pixels_i) - mean(green_pixels_j));
    diff_blue  = abs(mean(blue_pixels_i) - mean(blue_pixels_j));
    
    total_diff = sqrt(diff_red*diff_red + diff_green*diff_green + diff_blue*diff_blue);
    if(total_diff > 10)
        disocluded_regions(i) = 1;
    end
end

%% Visualization
disocluded_mask = zeros(s_x, s_y);
for i=1:s_x
    for j=1:s_y
        if(disocluded_regions(seg_i(i,j))==1)
            disocluded_mask(i,j) = 1;
        end
    end
end
figure;imshow(overlay_contour(I_i, disocluded_mask, [255 255 255]))
end

