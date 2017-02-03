function [partitions_vector_out] = create_clustered_partitions_sequence_intra(partitions_vector, regions_adjacency, partitions_adjacency, x)

[s_x, s_y, N_frames] = size(partitions_vector);

partitions_vector_out = zeros(s_x, s_y, N_frames);

for f=1:N_frames
    
    regions_partition = max(max(partitions_vector(:,:,f)));
    
    partition = ones(regions_partition, regions_partition);
    
    for i=1:length(regions_adjacency);
        
        partition_i = partitions_adjacency(1,i);
        partition_j = partitions_adjacency(2,i);
        
        if(partition_i==f && partition_j==f)
        
            region_i = regions_adjacency(1,i);
            region_j = regions_adjacency(2,i);

            if(x(i)==0)
                partition(region_i, region_j) = 0;
                partition(region_j, region_i) = 0;
            end
        
        end
    end

    labels = zeros(1, regions_partition);
    
    label=1;
    for ii=1:regions_partition
        if(labels(ii)==0)
            labels(ii)=label;
            labels = label_adjacents(ii, partition, labels);

            label = label+1;
        end
    end
    
    for i=1:s_x
        for j=1:s_y
            partitions_vector_out(i,j,f) = labels(partitions_vector(i,j,f));           
        end
    end
    
end

end