function [ edge_ids ] = select_edges_from_hierarchy( merging_sequence_i, RAG_i, regions_i, merging_sequence_j, RAG_j, regions_j, regions_adjacency)
%SELECT_EDGES_FROM_HIERARCHY Summary of this function goes here
%   Detailed explanation goes here
    edge_ids = [];
    
    for ii=1:size(merging_sequence_i,1)
        label1 = merging_sequence_i(ii,1);
        label2 = merging_sequence_i(ii,2);

        if label1 > regions_i || label2 > regions_i
            if label1 > regions_i
                label1 = get_subtree_leaves(merging_sequence_i,label1);
            end
            if label2 > regions_i
                label2 = get_subtree_leaves(merging_sequence_i,label2);
            end
            %find any region from label1 being adjacent to any region from
            %label2
            [row, col] = find(RAG_i(label1,label2)==1,1,'first');
            label1 = label1(row);
            label2 = label2(col);
        end
        
        if label1 > label2
            tmp = label1;
            label1 = label2;
            label2 = tmp;
        end
        
        %find label1-label2 pair in region_adjacency
        col_subset = find(regions_adjacency(1,:)==label1);
        edge_id = col_subset(find(regions_adjacency(2,col_subset)==label2));
        edge_ids = [edge_ids; edge_id];
    end
    
    for ii=1:size(merging_sequence_j,1)
        label1 = merging_sequence_j(ii,1);
        label2 = merging_sequence_j(ii,2);

        if label1 > regions_j || label2 > regions_j
            if label1 > regions_j
                label1 = get_subtree_leaves(merging_sequence_j,label1);
            end
            if label2 > regions_j
                label2 = get_subtree_leaves(merging_sequence_j,label2);
            end
            %find any region from label1 being adjacent to any region from
            %label2
            [row, col] = find(RAG_j(label1,label2)==1,1,'first');
            label1 = label1(row);
            label2 = label2(col);
        end
        
        if label1 > label2
            tmp = label1;
            label1 = label2;
            label2 = tmp;
        end
        
        label1 = label1 + regions_i;
        label2 = label2 + regions_i;
        %find label1-label2 pair in region_adjacency
        col_subset = find(regions_adjacency(1,:)==label1);
        edge_id = col_subset(find(regions_adjacency(2,col_subset)==label2));
        edge_ids = [edge_ids; edge_id];
    end
end

