function [partitions_vector_out, images_vector_out, last_label] = create_clustered_partitions_sequence_propagate(partitions_vector, partitions_vector_out, regions_adjacency, cumulate, x, labels_new1, labels_new2, last_label)
set(0,'RecursionLimit',2500)
[s_x, s_y, N_frames] = size(partitions_vector);

% total_regions = 0;
% for ii=1:N_frames
%     total_regions = total_regions + max(max(partitions_vector(:,:,ii)));
% end

total_regions = cumulate(end);

tmp = jet(double(total_regions+last_label));
color = uint8(255*tmp(randperm(total_regions+last_label,total_regions+last_label),:));

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

%% Assign labels of the first partition
for ii=1:length(labels_new1)
    labels(ii) = labels_new1(ii);
    
    labels = label_adjacents(ii, partition, labels);
end

for ii=1:length(labels_new2)
    if(labels(cumulate(2)+ii)==0)
        labels(cumulate(2)+ii) = labels_new2(ii);
        
        labels = label_adjacents(cumulate(2)+ii, partition, labels);
    end
end


label=last_label+1;
for ii=1:total_regions
    if(labels(ii)==0)
        labels(ii)=label;
        labels = label_adjacents(ii, partition, labels);
        
        label = label+1;
    end
end

last_label = label-1;

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