function [partition, merging_sequence] = ucm2bpt(ucm, max_num_leaves, old_format)
% [partition, merging_sequence] = ucm2bpt(ucm, max_num_leaves, old_format)
% ------------------------------------------------------------------------
%  ucm2bpt: 
%  Reads a UCM from file or directly as a matrix and converts it to a BPT
%   (partition and merging_sequence) where we accept non binary mergings,
%   that is, more than one merging for a given threshold.
%
%  Copyright (C)
%  Universitat Politecnica de Catalunya (UPC) - Barcelona - Spain
%
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  July 2012
% ------------------------------------------------------------------------

if(nargin<2)
    max_num_leaves = 0;
end
if(nargin<3)
    old_format = 1;
end

if ischar(ucm) % ucm refers to a file
    % Get UCM -> Must be saved as a variable named 'ucm2' or 'ucm.strength'
    load(ucm);
    if ~exist('ucm2', 'var')
        ucm2 = ucm.strength;
    end
elseif (size(ucm,1)>2 && size(ucm,2)>2) % It is a full ucm
    ucm2 = ucm;
    clear ucm;
else
    error('UCM type not accepted');
end
    
% Get leaves segmentation
labels = bwlabel(ucm2' == 0, 8); % Transposed for the scanning to be from
                                 %   left to right and from up to down
labels = labels';
segm = labels(2:2:end, 2:2:end);
num_leaves = max(segm(:));

% Get the merging thresholds
tmp = sort(unique(ucm2(:)));
ucm_thresholds = tmp(2:end);
num_mergings = length(ucm_thresholds);

% Store segmentation
partition = uint32(segm);

% ---------------------------


% Perform mergings
max_label = num_leaves;
curr_segm = segm;
total_mergings = 1;
total_merged_regs = 0;
is_binary = 1;
for ii=1:num_mergings
    curr_th=ucm_thresholds(ii);
    [x,y] = find(ucm2==curr_th);
    
    % Get the valid boundary positions
    b_idx = find(xor(mod(x,2),mod(y,2)));
    b_pos_x = x(b_idx);  b_pos_y = y(b_idx);
    
    % Separate them between up-down and left-right
    is_integ = logical(round((b_pos_x+1)/2)==(b_pos_x+1)/2);
    
    % Up-down positions
    label1_ud = curr_segm(sub2ind(size(curr_segm), (b_pos_x(is_integ)+1)/2, b_pos_y(is_integ)/2));
    label2_ud = curr_segm(sub2ind(size(curr_segm), (b_pos_x(is_integ)-1)/2, b_pos_y(is_integ)/2));        

    % Left-right positions
    label1_lr = curr_segm(sub2ind(size(curr_segm), b_pos_x(~is_integ)/2, (b_pos_y(~is_integ)+1)/2));
    label2_lr = curr_segm(sub2ind(size(curr_segm), b_pos_x(~is_integ)/2, (b_pos_y(~is_integ)-1)/2));

    labels = sort([label1_ud label2_ud;
                   label1_lr label2_lr],2);

    
    labels1 = labels(:,1);
    labels2 = labels(:,2);
    
    % Get 'unique' pairs
    K = 2*max_label+2; % 'infinity'
    ids = labels1 + K*labels2;
    ids = unique(ids);
    labels1 = mod(ids,K);
    labels2 = (ids-labels1)/K;
    
    % Detect non-binary fusions
    curr_merges = zeros(length(labels1), 2*length(labels1));
    curr_n_lab  = ones(length(labels1), 1);
    for kk=1:length(labels1)
        [xx1, yy1] = find(curr_merges==labels1(kk)); %#ok<NASGU>
        [xx2, yy2] = find(curr_merges==labels2(kk)); %#ok<NASGU>
        if isempty(xx1) && isempty(xx2)
            assert(curr_n_lab(kk)==1)
            curr_merges(kk,curr_n_lab(kk)) = labels1(kk);
            curr_merges(kk,curr_n_lab(kk)+1) = labels2(kk);
            curr_n_lab(kk) = 3;
        elseif isempty(xx1)
            is_binary = 0;
            curr_merges(xx2,curr_n_lab(xx2)) = labels1(kk);
            curr_n_lab(xx2) = curr_n_lab(xx2) + 1;
        elseif isempty(xx2)
            is_binary = 0;
            curr_merges(xx1,curr_n_lab(xx1)) = labels2(kk);
            curr_n_lab(xx1) = curr_n_lab(xx1) + 1;
        else
            is_binary = 0;
            if xx1~=xx2
                curr_merges(xx1,curr_n_lab(xx1):(curr_n_lab(xx1)+curr_n_lab(xx2)-2)) = ...
                    curr_merges(xx2,1:curr_n_lab(xx2)-1);
                curr_n_lab(xx1) = curr_n_lab(xx1)+curr_n_lab(xx2)-1;
                curr_merges(xx2,1:curr_n_lab(xx2)-1) = 0;
                curr_n_lab(xx2) = 1;
            end
        end
    end

    % Clean the result
    curr_merges = curr_merges(curr_n_lab~=1,:);
    curr_n_lab  = curr_n_lab(curr_n_lab~=1);
    
    % Perform mergings
    for kk=1:size(curr_merges,1)
        new_label = max_label + 1;
        max_label = new_label;

        for ll=1:curr_n_lab(kk)-1
            curr_segm(logical(curr_segm==curr_merges(kk,ll))) = new_label;            
        end

        merging_sequence{total_mergings}.sons = curr_merges(kk,1:curr_n_lab(kk)-1); %#ok<AGROW>
        merging_sequence{total_mergings}.father = new_label; %#ok<AGROW>

        total_mergings = total_mergings+1;
        total_merged_regs = total_merged_regs + curr_n_lab(kk)-2;
    end
end

if is_binary && old_format
    assert(total_mergings    == num_leaves);
    assert(total_merged_regs == num_leaves-1);
    
    % For backward compatibility, if tree is binary, store merging sequence in
    % the old way
    ms_tmp = zeros(num_leaves-1,3);
    for jj=1:num_leaves-1
        ms_tmp(jj,1) = merging_sequence{jj}.sons(1);
        ms_tmp(jj,2) = merging_sequence{jj}.sons(2);
        ms_tmp(jj,3) = merging_sequence{jj}.father;
    end
    merging_sequence = uint32(ms_tmp);
    
    % Prune and relabel if necessary
    if (max_num_leaves<num_leaves) && (max_num_leaves>0)
        [partition, merging_sequence] = prune_bpt(partition, merging_sequence, max_num_leaves);
        [partition, merging_sequence] = relabel_bpt(partition, merging_sequence);
    end
    
elseif ~is_binary && old_format
    error('Not implememented')
else
    assert(total_merged_regs == num_leaves-1);
    
    if (max_num_leaves<num_leaves) && (max_num_leaves>0)
        error('Not implememented')
    end
end





