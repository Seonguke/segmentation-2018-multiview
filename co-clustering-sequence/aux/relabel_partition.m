function [ relabeled_partition, labels_correspondence ] = relabel_partition( partition )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[sx, sy] = size(partition);
max_label = max(max(partition));

labels_new = zeros(1, max_label);

relabeled_partition = zeros(sx, sy);

label=1;
for i=1:sx
    for j=1:sy        
        if(labels_new(partition(i,j))==0)
            labels_new(partition(i,j)) = label;
            label = label + 1;
        end
        
        relabeled_partition(i,j) = labels_new(partition(i,j));
        
    end
end

labels_correspondence = zeros(1, label-1);
for i=1:sx
    for j=1:sy        
        labels_correspondence(relabeled_partition(i,j)) = partition(i,j);
    end
end


end

