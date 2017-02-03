function [ disoccluded_mask ] = disocclusion_detection( I_i, I_j, seg_i, seg_j, OF_regions_x, OF_regions_y )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[s_x, s_y] = size(seg_j);
predicted_image = zeros(s_x, s_y, 3);
counter_estimated_pixels = zeros(s_x, s_y);

%% Create predicted image using OF
for i=1:s_x
    for j=1:s_y        
        delta_x = round(OF_regions_x(seg_i(i,j)));
        delta_y = round(OF_regions_y(seg_i(i,j)));
        
        %% Accumulate pixel values
        if((i+delta_y)>0 && (i+delta_y)<=s_x && (j+delta_x)>0 && (j+delta_x)<=s_y)
            predicted_image(i+delta_y,j+delta_x, 1) = predicted_image(i+delta_y,j+delta_x, 1) + I_i(i,j,1);  
            predicted_image(i+delta_y,j+delta_x, 2) = predicted_image(i+delta_y,j+delta_x, 2) + I_i(i,j,2);  
            predicted_image(i+delta_y,j+delta_x, 3) = predicted_image(i+delta_y,j+delta_x, 3) + I_i(i,j,3);  
        
            counter_estimated_pixels(i+delta_y,j+delta_x) = counter_estimated_pixels(i+delta_y,j+delta_x) + 1;
        end      
    end
end
%% Compute de average color for each pixel
for i=1:s_x
    for j=1:s_y        
        if(counter_estimated_pixels(i,j)>0)
            predicted_image(i,j,1) = predicted_image(i,j,1)/counter_estimated_pixels(i,j);
            predicted_image(i,j,2) = predicted_image(i,j,2)/counter_estimated_pixels(i,j);
            predicted_image(i,j,3) = predicted_image(i,j,3)/counter_estimated_pixels(i,j);  
        end
    end
end
% figure;imshow(uint8(predicted_image))


n_regions = max(max(seg_j));
mean_color_regions = zeros(3,n_regions);
mean_color_prediction = zeros(3,n_regions);
area_regions = zeros(1, n_regions);
%% Compute the mean color of regions using the original and the predicted images
for i=1:s_x
    for j=1:s_y        
        mean_color_regions(1,seg_j(i,j)) = mean_color_regions(1,seg_j(i,j)) + double(I_j(i,j,1));
        mean_color_regions(2,seg_j(i,j)) = mean_color_regions(2,seg_j(i,j)) + double(I_j(i,j,2));
        mean_color_regions(3,seg_j(i,j)) = mean_color_regions(3,seg_j(i,j)) + double(I_j(i,j,3));
        
        mean_color_prediction(1,seg_j(i,j)) = mean_color_prediction(1,seg_j(i,j)) + predicted_image(i,j,1);
        mean_color_prediction(2,seg_j(i,j)) = mean_color_prediction(2,seg_j(i,j)) + predicted_image(i,j,2);
        mean_color_prediction(3,seg_j(i,j)) = mean_color_prediction(3,seg_j(i,j)) + predicted_image(i,j,3);
        
        area_regions(seg_j(i,j)) = area_regions(seg_j(i,j)) + 1;
    end
end

mean_color_regions(1,:) = round(mean_color_regions(1,:)./area_regions);
mean_color_regions(2,:) = round(mean_color_regions(2,:)./area_regions);
mean_color_regions(3,:) = round(mean_color_regions(3,:)./area_regions);  

mean_color_prediction(1,:) = round(mean_color_prediction(1,:)./area_regions);
mean_color_prediction(2,:) = round(mean_color_prediction(2,:)./area_regions);
mean_color_prediction(3,:) = round(mean_color_prediction(3,:)./area_regions); 

%% Visualize regions and mean colors
mean_regions_image = zeros(s_x, s_y, 3);
mean_regions_prediction = zeros(s_x, s_y, 3);
for i=1:s_x
    for j=1:s_y        
        mean_regions_image(i,j,1) = mean_color_regions(1,seg_j(i,j));
        mean_regions_image(i,j,2) = mean_color_regions(2,seg_j(i,j));
        mean_regions_image(i,j,3) = mean_color_regions(3,seg_j(i,j));
        
        mean_regions_prediction(i,j,1) = mean_color_prediction(1,seg_j(i,j));
        mean_regions_prediction(i,j,2) = mean_color_prediction(2,seg_j(i,j));
        mean_regions_prediction(i,j,3) = mean_color_prediction(3,seg_j(i,j));
    end
end

% figure;imshow(uint8(mean_regions_image))
% figure;imshow(uint8(mean_regions_prediction))

%% Compute the difference between the mean colors of regions 
%% evaluated over the original and the predicted images
th = 80;
diff = abs(mean_regions_image - mean_regions_prediction);
distance = sqrt(diff(:,:,1).*diff(:,:,1) + diff(:,:,2).*diff(:,:,2) + diff(:,:,3).*diff(:,:,3));
% figure;imshow(distance > th)

disoccluded_mask = distance > th;
end

