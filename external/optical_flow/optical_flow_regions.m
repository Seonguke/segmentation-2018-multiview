function [ OF_regions_x, OF_regions_y, flow ] = optical_flow_regions( I_i, I_j, seg_i, seg_j )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% [ I_i, seg_i, vectors_i ] = read_measurement( sequence_path, sequence_name, extension, seq_info.frame -1, th_ucm );
% [ I_j, seg_j, vectors_j ] = read_measurement( sequence_path, sequence_name, extension, seq_info.frame, th_ucm );
% 
% show(I_i, overlay_contour(label2rgb(seg_i), seg_i, [0 0 0]), debug);
% show(I_j, overlay_contour(label2rgb(seg_j), seg_j, [0 0 0]), debug);

    
%%%% Compute LDOF of the images %%%%
flow = mex_LDOF(double(I_i),double(I_j));
%%%%%%%%%%%%%%%%%%%

%%%% Average of the regions OF %%%%
n_regions = max(max(seg_i));
OF_regions_x = zeros(1,n_regions);
OF_regions_y = zeros(1,n_regions);
area_regions = zeros(1,n_regions);

[s_x, s_y] = size(seg_i);
for i=1:s_x
    for j=1:s_y
        OF_regions_x(seg_i(i,j)) = OF_regions_x(seg_i(i,j)) + flow(i,j,1);
        OF_regions_y(seg_i(i,j)) = OF_regions_y(seg_i(i,j)) + flow(i,j,2);        
        
        area_regions(seg_i(i,j)) = area_regions(seg_i(i,j)) + 1; 
    end
end

OF_regions_x = OF_regions_x./area_regions;
OF_regions_y = OF_regions_y./area_regions;

%% Visualize Optical Flow
% centroids  = regionprops(seg_i, 'centroid');
% u = zeros(1, length(centroids));
% v = zeros(1, length(centroids));
% 
% for i=1:length(centroids)
% u(i) = centroids(i).Centroid(1);
% v(i) = centroids(i).Centroid(2);
% end
% 
% figure;imshow(I_i);hold on;quiver(u, v, OF_regions_x, OF_regions_y);
% 
% u = 1:3:s_x;
% v = 1:3:s_y;
% [a,b]=meshgrid(v,u);
% 
% OF_regions_x = zeros(1,length(u));
% OF_regions_y = OF_regions_x;
% 
% for i=1:length(u)
%     OF_regions_x(i) = flow(u(i),v(i),1);
%     OF_regions_x(i) = flow(u(i),v(i),2);
% end
% 
% % figure;imshow(I_i);hold on;quiver(u, v, flow([u' v'],1), flow(u,v,2));
% figure;imshow(I_i);hold on;quiver(a, b, flow(u,v,1), flow(u,v,2));
% 
% figure;imshow(I_j);
%%

