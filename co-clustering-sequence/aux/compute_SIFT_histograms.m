function [ SIFT ] = compute_SIFT_histograms( I, labeled_elements )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
elements = max(max(labeled_elements));
%scale = 8/3;
scale = 16/3;
%scale = 24/3;
%scale = 32/3;

%scale = 4/3;

elements_positions = find_elements_position(labeled_elements, elements);
elements_positions = elements_positions/2;
%SIFT = zeros(128,elements);
% for i=1:elements
%     x1 = elements_positions(i,1);
%     y1 = elements_positions(i,2);
%     d = compute_SIFT_elem( I, x1, y1, scale );
%     SIFT(:,i) = d;
% end
frames = [elements_positions, scale*ones(size(elements_positions,1),1), zeros(size(elements_positions,1),1)];
[f1,SIFT]=vl_sift(im2single(rgb2gray(I)), 'frames', frames');
SIFT = double(SIFT);
normSIFT = repmat(sum(SIFT,1),128,1);
normSIFT(:,find(sum(SIFT,1)==0))=1;
SIFT = SIFT./normSIFT;

end