function [ucm2, leaves_partition, merging_sequence] = im2ucm(img_file, ucm_file, bpt_file, leaves_part_file) %#ok<STOUT>

if nargin==3
    error('To write a BPT to file, both "bpt_file" and "leaves_part_file" must be provided')
end

% Call globalPb using evalc to avoid all the messages that are written
evalc(['gPb_orient = globalPb(''' img_file ''',''' ucm_file ''')']);
evalc('ucm2 = contours2ucm(gPb_orient, ''doubleSize'')');

if nargin>1
    save(ucm_file,'ucm2');
end

if nargout>2 || nargin>3
    [leaves_partition, merging_sequence] = ucm2bpt(ucm2);
end

if nargin>3
    write_bpt(bpt_file,leaves_part_file,img_file,leaves_partition,merging_sequence,1)
end




