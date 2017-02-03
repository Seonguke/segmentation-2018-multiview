function [rel_part, rel_mer_seq] = relabel_bpt(partition, merging_sequence)
% [rel_part, rel_mer_seq] = relabel_bpt(partition, merging_sequence)
% ------------------------------------------------------------------------
%  relabel_bpt: 
%  Relabels a BPT, partition from 1 to N, mergings from N+1 to 2N-1
%
%  Copyright (C)
%  Universitat Politecnica de Catalunya (UPC) - Barcelona - Spain
%
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2011
% ------------------------------------------------------------------------

% Relabel partition
rel_part = uint32(relabel(partition));

% Create LUT between relabeled and original
K = max( max(partition(:)), max(rel_part(:)))+1;
idx_matrix = partition+K*rel_part;
idx = unique(idx_matrix);
old_labels = mod(idx,K);
new_labels = (idx-old_labels)/K;

% Relabel merging sequence
rel_mer_seq = zeros(size(merging_sequence));
nm = size(merging_sequence,1);
nr = length(old_labels);
for ii=1:nm
    rel_mer_seq(ii,1) = new_labels(find(old_labels==merging_sequence(ii,1), 1, 'first'));
    rel_mer_seq(ii,2) = new_labels(find(old_labels==merging_sequence(ii,2), 1, 'first'));
    rel_mer_seq(ii,3) = nr+ii;
    
    old_labels = [old_labels; merging_sequence(ii,3)]; %#ok<AGROW>
    new_labels = [new_labels; nr+ii]; %#ok<AGROW>
end

% to uint32
rel_part = uint32(rel_part);
rel_mer_seq = uint32(rel_mer_seq);
