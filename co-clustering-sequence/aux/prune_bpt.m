function [pru_part, pru_mer_seq] = prune_bpt(partition, merging_sequence, num_leaves)
% [pru_part, pru_mer_seq] = prune_bpt(partition, merging_sequence, num_leaves)
% ------------------------------------------------------------------------
%  prune_bpt: 
%  Prunes a BPT at a given number of leaves
%
%  Copyright (C)
%  Universitat Politecnica de Catalunya (UPC) - Barcelona - Spain
%
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2011
% ------------------------------------------------------------------------

nr = length(unique(partition));
num_mergs_to_do = nr-num_leaves;
assert(num_mergs_to_do>=0)
assert(num_mergs_to_do<nr)

% Get leaves partition
pru_part = partition;
for ii=1:num_mergs_to_do
    % Perform merging
    pru_part(logical(pru_part==merging_sequence(ii,1))) = merging_sequence(ii,3);
    pru_part(logical(pru_part==merging_sequence(ii,2))) = merging_sequence(ii,3);
end

% Cut merging_sequence
pru_mer_seq = merging_sequence(num_mergs_to_do+1:end,:);
pru_part = uint32(pru_part);
