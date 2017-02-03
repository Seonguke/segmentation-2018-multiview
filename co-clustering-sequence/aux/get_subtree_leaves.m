function sel_regs = get_subtree_leaves( merging_sequence, selected_region)
%
% ------------------------------------------------------------------------
%  Copyright (C)
%  Universitat Politecnica de Catalunya (UPC) - Barcelona - Spain
%
%  Jordi Pont-Tuset <jordi.pont@upc.edu>
%  June 2011
% ------------------------------------------------------------------------
sel_regs = selected_region;
for ii=size(merging_sequence,1):-1:1
    where_is = find(merging_sequence(ii,3)==sel_regs,1);
    if ~isempty(where_is)
        sel_regs(where_is) = [];
        sel_regs = [sel_regs  merging_sequence(ii,1:2)]; %#ok<AGROW>
    end
end
sel_regs = sort(sel_regs);
