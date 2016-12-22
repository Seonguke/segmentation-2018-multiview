function [partitions_vector_out, images_vector_out] = create_clustered_partitions_sequence_global_opt(partitions_vector, regions_adjacency, cumulate, x)
set(0,'RecursionLimit',2500)
[s_x, s_y, N_frames] = size(partitions_vector);

total_regions = cumulate(end);

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

partitions_vector_out = zeros(s_x, s_y, N_frames);
images_vector_out     = zeros(s_x, s_y, 3, N_frames);

for i=1:s_x
    for j=1:s_y
        
        for k=1:N_frames
            partitions_vector_out(i,j,k) = labels(cumulate(k)+partitions_vector(i,j,k));
            
            images_vector_out(i,j,1,k) = color(labels(cumulate(k)+partitions_vector(i,j,k)),1);
            images_vector_out(i,j,2,k) = color(labels(cumulate(k)+partitions_vector(i,j,k)),2);
            images_vector_out(i,j,3,k) = color(labels(cumulate(k)+partitions_vector(i,j,k)),3);
        
        end

    end
end

end