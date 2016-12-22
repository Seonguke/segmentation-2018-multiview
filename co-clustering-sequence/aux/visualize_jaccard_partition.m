function [ partition ] = visualize_jaccard_partition( x, seg )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[s_x, s_y] = size(seg);
partition = ones(s_x,s_y);

region_count = 2;

for ii = 1:length(x)
    if(x(ii)==1)
            partition = region_count*double(seg==ii)+partition.*double(seg~=ii);
            region_count = region_count + 1;
    end
end
end