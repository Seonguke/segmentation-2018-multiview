function [partition_i, partition_j] = create_clustered_partitions_original2_intra(seg_i, seg_j, regions_adjacency, x)
regions_i = max(max(seg_i));
regions_j = max(max(seg_j));

partition_i_out = ones(regions_i, regions_i);
partition_j_out = ones(regions_j, regions_j);

for i=1:length(regions_adjacency);
    region_i = regions_adjacency(1,i);
    region_j = regions_adjacency(2,i);
    
    if(region_i<=regions_i && region_j<=regions_i)
    
        if(x(i)==0)
            partition_i_out(region_i, region_j) = 0;
            partition_i_out(region_j, region_i) = 0;
        end   
        
    elseif(region_i>regions_i && region_j>regions_i)
        
        if(x(i)==0)
            partition_j_out(region_i-regions_i, region_j-regions_i) = 0;
            partition_j_out(region_j-regions_i, region_i-regions_i) = 0;
        end
        
    end
end


labels_i = zeros(1, regions_i);

label=1;
for ii=1:regions_i
    if(labels_i(ii)==0)
        labels_i(ii)=label;
        labels_i = label_adjacents(ii, partition_i_out, labels_i);
        
        label = label+1;
    end
end

labels_j = zeros(1, regions_j);label=1;
for ii=1:regions_j
    if(labels_j(ii)==0)
        labels_j(ii)=label;
        labels_j = label_adjacents(ii, partition_j_out, labels_j);
        
        label = label+1;
    end
end
        
s1 = size(seg_i,1);
s2 = size(seg_i,2);

partition_i = zeros(s1,s2);
partition_j = zeros(s1,s2);


for i=1:s1
    for j=1:s2
        
        partition_i(i,j) = labels_i(seg_i(i,j));
        partition_j(i,j) = labels_j(seg_j(i,j));
        
    end
end

end