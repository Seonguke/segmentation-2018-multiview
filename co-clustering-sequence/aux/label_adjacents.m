function labels = label_adjacents( region, partition, labels )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
total_regions =length(labels);

new_adjacents = zeros(1, total_regions);

for ii=1:total_regions
    if(partition(region,ii)==0 && labels(ii)==0)
        labels(ii) = labels(region);     
        new_adjacents(ii) = 1;
    end
end

for ii=1:total_regions
    if(new_adjacents(ii)==1)
        labels = label_adjacents(ii, partition, labels);
    end
end
    
end