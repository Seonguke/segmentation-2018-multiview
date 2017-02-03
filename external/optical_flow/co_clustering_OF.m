function [ x, partition1, partition2, Q, W_count ] = co_clustering_OF( I_i, I_j, seg_i, seg_j, vectors_i, vectors_j, seq_info )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
% show(label2rgb(seg_i),label2rgb(seg_j), debug);

%%%% Obtain labelled elements %%%%
[ labeled_elements_i, gridbmap_i, idx_neighbors_i, angles_i ] = obtain_contour_elements( seg_i, vectors_i );
[ labeled_elements_j, gridbmap_j, idx_neighbors_j, angles_j ] = obtain_contour_elements( seg_j, vectors_j );
%%%%%%%%%%%%%%%%%%%%%%%

%%%% Compute Bi matrix %%%%
[ B_i, gPb_angles_i ] = computeBmatrix(I_i, seg_i, angles_i, labeled_elements_i, idx_neighbors_i.matrix_max, idx_neighbors_i.matrix_min);
% image_angles_i = zeros(size(gPb_angles_i,1), size(gPb_angles_i,2));
% colors = floor(linspace(50,255,8));
% for i=1:size(gPb_angles_i,1)
%     for j=1:size(gPb_angles_i,2)
%     if(gPb_angles_i(i,j)==22.5)
%         image_angles_i(i,j)=colors(1);
%     elseif(gPb_angles_i(i,j)==45)
%         image_angles_i(i,j)=colors(2);
%     elseif(gPb_angles_i(i,j)==45+22.5)
%         image_angles_i(i,j)=colors(3);
%     elseif(gPb_angles_i(i,j)==90)
%         image_angles_i(i,j)=colors(4);
%     elseif(gPb_angles_i(i,j)==90+22.5)
%         image_angles_i(i,j)=colors(5);
%     elseif(gPb_angles_i(i,j)==90+45)
%         image_angles_i(i,j)=colors(6);
%     elseif(gPb_angles_i(i,j)==90+45+22.5)
%         image_angles_i(i,j)=colors(7);
%     elseif(gPb_angles_i(i,j)==180)
%         image_angles_i(i,j)=colors(8);
%     end
%     end
% end
% CORRECTO
% Elapsed time 2.2 seconds
% Se puede pasar a mex
%%%%%%%%%%%%%%%%%%%%%%%

%%%% Load Bj matrix and gPb_anglesj %%%%
% [ B_j, gPb_angles_j ] = computeBmatrix(I_j, seg_j, angles_j, labeled_elements_j, idx_neighbors_j.matrix_max, idx_neighbors_j.matrix_min);
[ B_j, gPb_angles_j ] = load_B_angles( seq_info );
% image_angles_j = zeros(size(gPb_angles_i,1), size(gPb_angles_i,2));
% colors = floor(linspace(50,255,8));
% for i=1:size(gPb_angles_j,1)
%     for j=1:size(gPb_angles_j,2)
%     if(gPb_angles_j(i,j)==22.5)
%         image_angles_j(i,j)=colors(1);
%     elseif(gPb_angles_j(i,j)==45)
%         image_angles_j(i,j)=colors(2);
%     elseif(gPb_angles_j(i,j)==45+22.5)
%         image_angles_j(i,j)=colors(3);
%     elseif(gPb_angles_j(i,j)==90)
%         image_angles_j(i,j)=colors(4);
%     elseif(gPb_angles_j(i,j)==90+22.5)
%         image_angles_j(i,j)=colors(5);
%     elseif(gPb_angles_j(i,j)==90+45)
%         image_angles_j(i,j)=colors(6);
%     elseif(gPb_angles_j(i,j)==90+45+22.5)
%         image_angles_j(i,j)=colors(7);
%     elseif(gPb_angles_j(i,j)==180)
%         image_angles_j(i,j)=colors(8);
%     end
%     end
% end
% Elapsed time 0.11 seconds
%%%%%%%%%%%%%%%%%%%%%%%

%%%% Compute descriptors %%%%

H_i = compute_histograms( I_i, gridbmap_i, labeled_elements_i );
H_j = compute_histograms( I_j, gridbmap_j, labeled_elements_j );
% Elapsed time 5.93 seconds
% Depende de los elementos de contorno

HOG_i = compute_HOG_histograms( gPb_angles_i, gridbmap_i, labeled_elements_i );
HOG_j = compute_HOG_histograms( gPb_angles_j, gridbmap_j, labeled_elements_j );
% Elapsed time 0.56 seconds
% Depende de los elementos de contorno

d_i = [H_i; HOG_i]';
d_j = [H_j; HOG_j]';
%%%%%%%%%%%%%%%%%%%%


%%%% Compute Q matrix %%%%
Q11 = compute_Q_intra(I_i, seg_i, idx_neighbors_i.matrix_max, idx_neighbors_i.matrix_min);
Q22 = compute_Q_intra(I_j, seg_j, idx_neighbors_j.matrix_max, idx_neighbors_j.matrix_min);
% Elapsed time 0.56 seconds

%5,3.5
sigma_H   = 60;
sigma_HOG = 0;
%13.5;
Inv_Sigma = create_sigma_matrix(H_i, HOG_i, sigma_H, sigma_HOG);

[W, W_count] = compute_Q_inter(d_i, d_j, labeled_elements_i, labeled_elements_j, Inv_Sigma);
% Elapsed time 0.39 seconds

Q12 = B_i'*W*B_j;
Q21 = Q12';

Q = [Q11 Q12; Q21 3*Q22]; 
W_count =  double(B_i~=0)'*W_count*double(B_j~=0);
%%%%%%%%%%%%%%%%%%%%

%%%% Compute Adjacency Graph %%%%
[ RAG_i, RAG_j ] = compute_adjacency_graph( seg_i, seg_j );
% Elapsed time 0.19 seconds

[ OF_regions_x, OF_regions_y, flow ] = optical_flow_regions( I_i, I_j, seg_i, seg_j );
% Elapsed time 11.34 seconds

RAG_ij = compute_adjacency_OF( seg_i, seg_j, OF_regions_x, OF_regions_y );
% Elapsed time 0.009 seconds

RAG = [RAG_i RAG_ij;
       RAG_ij' RAG_j];
%%%%%%%%%%%%%%%%%%%%

%%% Variables creation %%%%
regions_i = max(max(idx_neighbors_i.matrix_max));
regions_j = max(max(idx_neighbors_j.matrix_max));
total_regions = regions_i + regions_j;
variables =(total_regions)*(total_regions-1)/2;
%%%%%%%%%%%%%%%%%%%%

%%%% Create Constraints Matrix %%%%
% [A, b, Aeq, beq] = constraints(regions_i, regions_j, RAG_i, RAG_j, variables);
[A, b, Aeq, beq] = constraints_OF(regions_i, regions_j, RAG, variables);
% Elapsed time 0.47 seconds
%%%%%%%%%%%%%%%%%%%%

%%%% ub & lb %%%%
ub = ones(1, variables);
lb = zeros(1, variables);
%%%%%%%%%%%%%%%%%%%%

% f=[];
f = zeros(1, variables);
count = 1;
for i=1:total_regions
    for j=i+1:total_regions
%         f = [f Q(i,j)+Q(j,i)];
        f(count) = Q(i,j)+Q(j,i);
        count = count + 1;
    end
end
% Elapsed time 0.01 seconds

addpath('/usr/local/opt/CPLEX_Studio124/cplex/matlab/')
% x = cplexlp(f',A,b',Aeq,beq',lb',ub');

x = cplexbilp(f',A,b',Aeq,beq');
% Elapsed time 0.41 seconds

[partition1, partition2] = create_clustered_partitions_OF2(I_i, I_j, seg_i, seg_j, RAG, x);
%%%%%%%%%%%%%%%%%%%%%%%
end