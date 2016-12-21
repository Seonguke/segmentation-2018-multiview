function [ adjacency_matrix ] = compute_adjacency_OF( seg_i, seg_j, OF_regions_x, OF_regions_y, object_mask, flow )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Variables declaration
regions_i = max(max(seg_i));
regions_j = max(max(seg_j));

adjacency_matrix_with_OF = zeros(regions_i, regions_j);
adjacency_matrix_without_OF = zeros(regions_i, regions_j);

area_regions_j = zeros(1,regions_j);
area_regions_i = zeros(1,regions_i);

[s_x, s_y] = size(seg_i);

%% Compute adjacency and areas
for i=1:s_x
    for j=1:s_y        
        %% Compute adjacent area with OF
        delta_x = round(OF_regions_x(seg_i(i,j)));
        delta_y = round(OF_regions_y(seg_i(i,j)));
        
        if((i+delta_y)>0 && (i+delta_y)<=s_x && (j+delta_x)>0 && (j+delta_x)<=s_y)
            adjacency_matrix_with_OF(seg_i(i,j),seg_j(i+delta_y,j+delta_x)) = adjacency_matrix_with_OF(seg_i(i,j),seg_j(i+delta_y,j+delta_x)) + 1;           
        end      
        %% Compute adjacent area without OF
        adjacency_matrix_without_OF(seg_i(i,j),seg_j(i,j)) = adjacency_matrix_without_OF(seg_i(i,j),seg_j(i,j)) + 1;
        
        %% Compute region areas
        area_regions_i(seg_i(i,j)) = area_regions_i(seg_i(i,j)) + 1;
        area_regions_j(seg_j(i,j)) = area_regions_j(seg_j(i,j)) + 1; 
    end
end
%% Get regions from mask
regions_i = max(max(seg_i));
regions_j = max(max(seg_j));

regions_mask = zeros(1, regions_i);
for i=1:s_x
    for j=1:s_y
        if(object_mask(i,j))
            regions_mask(seg_i(i,j)) = 1;
        end
    end
end
%% Compute regions adjacency with OF
for j=1:regions_j
    if(j==231)
        1
    end
    total_area_adjacency = 0;
    adjacents_from_object = zeros(regions_i,1);
    %% Minimum area of object regions included in region from seg_j
    for i=1:regions_i
        if(regions_mask(i) == 1)
            if(adjacency_matrix_with_OF(i,j)/area_regions_i(i) > 0.02)
                total_area_adjacency = total_area_adjacency + adjacency_matrix_with_OF(i,j);
                adjacents_from_object(i) = 1;
            end
        end
    end
    %% Minimum area of regions from seg_j covered by object regions
    if(total_area_adjacency/area_regions_j(j) > 0.1)
        adjacency_matrix_with_OF(:,j) = adjacents_from_object;
    else
        adjacency_matrix_with_OF(:,j) = (adjacency_matrix_with_OF(:,j)>0).*(~regions_mask)';
    end
end
    
%% Compute regions adjacency without OF
for j=1:regions_j
%     adjacency_matrix_with_OF(:,i) = adjacency_matrix_with_OF(:,i)/area_regions(i);
    adjacency_matrix_without_OF(:,j) = adjacency_matrix_without_OF(:,j)/area_regions_j(j);
end
%% Only consider as adjacent those regions with a minimum percentage
th_s = 0.2;
% adjacency_matrix_with_OF = adjacency_matrix_with_OF > th_s;
adjacency_matrix_without_OF = adjacency_matrix_without_OF > th_s;

adjacency_matrix = adjacency_matrix_with_OF | adjacency_matrix_without_OF;

%% Visualize Optical Flow
% translated_regions = zeros(s_x, s_y);
% for i=1:s_x
%     for j=1:s_y
%         if((i + delta_x)>0 &&  (j + delta_y)>0 && (i + delta_x)<= s_x &&  (j + delta_y)<= s_y)
%             translated_regions(i + delta_x, j + delta_y) = seg_i(i,j);
%         end
%     end
% end
end

