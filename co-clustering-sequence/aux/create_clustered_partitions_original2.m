function [partition_i, partition_j, image_i, image_j] = create_clustered_partitions_original2(seg_i, seg_j, regions_adjacency, x)
regions_i = max(max(seg_i));
regions_j = max(max(seg_j));

total_regions = regions_i + regions_j;

tmp = jet(double(total_regions));
color = uint8(255*tmp(randperm(total_regions,total_regions),:));

partition = ones(total_regions, total_regions);
for i=1:length(regions_adjacency);
    region_i = regions_adjacency(1,i);
    region_j = regions_adjacency(2,i);
    
    if(x(i)==0)
        partition(region_i, region_j) = 0;
        partition(region_j, region_i) = 0;
    end
end


labels = zeros(1, total_regions);

label=1;
for ii=1:total_regions
    if(labels(ii)==0)
        labels(ii)=label;
        labels = label_adjacents(ii, partition, labels);
        
        label = label+1;
    end
end
        
s1 = size(seg_i,1);
s2 = size(seg_i,2);

partition_i = zeros(s1,s2);
partition_j = zeros(s1,s2);

image_i = zeros(s1,s2,3);
image_j = zeros(s1,s2,3);

for i=1:s1
    for j=1:s2
        
        partition_i(i,j) = labels(seg_i(i,j));
        partition_j(i,j) = labels(regions_i+seg_j(i,j));
        
        image_i(i,j,1) = color(labels(seg_i(i,j)),1);
        image_i(i,j,2) = color(labels(seg_i(i,j)),2);
        image_i(i,j,3) = color(labels(seg_i(i,j)),3);
        
        image_j(i,j,1) = color(labels(regions_i+seg_j(i,j)),1);
        image_j(i,j,2) = color(labels(regions_i+seg_j(i,j)),2);
        image_j(i,j,3) = color(labels(regions_i+seg_j(i,j)),3);

    end
end

end