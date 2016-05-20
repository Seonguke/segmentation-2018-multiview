function [ nodes_leaves, variables_intra ] = compute_nodes_leaves( partition, merging_sequence )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n_fusions = length(merging_sequence);
regions_partition = max(max(partition));
nodes = regions_partition + n_fusions -1;

nodes_leaves = zeros(nodes, nodes);
variables_intra = zeros(2, n_fusions-1);

for i=1:regions_partition
    nodes_leaves(i,i) = 1;
end

for i=1:(n_fusions-1)
    reg_A = merging_sequence(i,1);
    reg_B = merging_sequence(i,2);
    reg_C = merging_sequence(i,3);
    
    nodes_leaves(:,reg_C) = nodes_leaves(:,reg_A) | nodes_leaves(:,reg_B);
    
    nodes_leaves(reg_C,reg_C) = 1;
    
    variables_intra(1,i) = reg_A;
    variables_intra(2,i) = reg_B;
end


end

