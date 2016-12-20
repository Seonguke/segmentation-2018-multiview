function [partition_i_inter, partition_j_inter, image_i_inter, image_j_inter, partition_i_intra, partition_j_intra] = ...
    co_clustering_global_hierarchy_optical_flow_num_clusters( sequence_name, idx_j, I_i, I_j, seg_i, seg_j, ucm_i, ucm_j, merging_sequence_i, merging_sequence_j, vectors_i, vectors_j, window_similarities, w_color, w_hog, w_sift, hierarchical_constraint_opt, num_clusters  )
% [ partition_i, partition_j, x, Q ] = co_clustering_original( I_i, I_j, seg_i, seg_j, vectors_i, vectors_j )
% ------------------------------------------------------------------------
% Performs an unsupervised shape-based co-clustering of two segmentations.
%
% INPUT
%	I_i, I_j                    Input images.
%	seg_i, seg_j                Input image segmentations.
%   vectors_i, vectors_j        Contour vectors resulting from the function 
%                               create_finest_partition (performed during 
%                               the ucm creation).
%
% OUTPUT
%	partition_i, partition_j	Co-clustered partitions.
%   Q                           Affinity matrix between regions.
%   x                           Vector of distances between regions.
%                               0 indicates neightbor regions, 1 indicates 
%                               that regions belong to different clusters.
%
% ------------------------------------------------------------------------
%  Copyright (C)
%  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
%
%  David Varas <david.varas@upc.edu>
%  May 2014
% ------------------------------------------------------------------------

%% Obtain labeled elements 
% Label contour elements in order to keep 
% the coherence between functions
% fprintf(' Label contour elements...  \n')
%tic;
[ labeled_elements_i, gridbmap_i, idx_neighbors_i, angles_i ] = obtain_contour_elements(double(seg_i), vectors_i );
[ labeled_elements_j, gridbmap_j, idx_neighbors_j, angles_j ] = obtain_contour_elements(double(seg_j), vectors_j );
%toc;
% fprintf(' ...ended! \n')
%%%%%%%%%%%%%%%%%%%%%%%

%% Compute basic information 
% Number of regions in the partitions 
% and number of variables in the process
fprintf(' Compute basic information... \n')
%tic;
regions_i = max(max(seg_i));
regions_j = max(max(seg_j));
%total_regions = regions_i+regions_j;

%toc;
fprintf(' ...ended! \n')
%%%%%%%%%%%%%%%%%%%%%%%

%% Compute adjacency between partitions (agirbau)
% This adjacency is used to reduce the 
% number of variables in the process
% fprintf(' Compute Adjacency... \n')
%tic;
RAG_i = RAG_Intra_mex(idx_neighbors_i,double(seg_i));
RAG_j = RAG_Intra_mex(idx_neighbors_j,double(seg_j));

load(['/work/cventura/segmentation_propagation/optical_flow/' sequence_name '/frame_OF_' sprintf('%05d',idx_j) '.mat']);
MV_x = flow(:,:,1);
MV_y = flow(:,:,2);
RAG_ij = RAG_Inter_motion_estimation_optical_flow_mex(double(seg_i), double(seg_j), MV_x, MV_y);


RAG =[RAG_i RAG_ij;
      RAG_ij' RAG_j];

variables_opt = sum(sum(RAG==1))/2;
%toc;
% fprintf(' ...ended! \n')
%%%%%%%%%%%%%%%%%%%%%%%

%% Variables selection 
% Variables are stored in regions_adjacency
% The variable that encodes the adjacency 
% between regions i and j is stored in regions_variables(i,j)
% fprintf(' Variables selection... \n')
%tic;
[regions_adjacency, regions_variables] = variables_reduction(RAG, variables_opt);
%toc;
% fprintf(' ...ended! \n')
%%%%%%%%%%%%%%%%%%%%%%%

%% Compute B matrices (agirbau)
% B matrices contain the vectors normal to 
% the contour at each contour element pointing 
% outside the region
% fprintf(' Compute B matrices...  \n')
% tic;
[B_i, gPb_angles_i] = computeBmatrix_sparse_mex(double(seg_i),angles_i,labeled_elements_i);
[B_j, gPb_angles_j] = computeBmatrix_sparse_mex(double(seg_j),angles_j,labeled_elements_j);
% [B_i, gPb_angles_i] = computeBmatrix_mex(double(seg_i),angles_i,labeled_elements_i);
% [B_j, gPb_angles_j] = computeBmatrix_mex(double(seg_j),angles_j,labeled_elements_j);


%toc;
% fprintf(' ...ended! \n')
%%%%%%%%%%%%%%%%%%%%%%%

%% Compute descriptors (agirbau)
% Compute descritors for each contour element
% fprintf(' Compute descriptors...  \n')
%tic;
[d_i,n_bins_RGB_i,n_bins_grad_i] = calculate_descriptors_mex(double(I_i),gPb_angles_i,labeled_elements_i, double(idx_neighbors_i.matrix_min), double(idx_neighbors_i.matrix_max));
[d_j,n_bins_RGB_j,n_bins_grad_j] = calculate_descriptors_mex(double(I_j),gPb_angles_j,labeled_elements_j, double(idx_neighbors_i.matrix_min), double(idx_neighbors_i.matrix_max));
%toc;
% fprintf(' ...ended! \n')
%%%%%%%%%%%%%%%%%%%%

[ d_i2 ] = compute_HOG_descriptors( I_i, labeled_elements_i );
[ d_j2 ] = compute_HOG_descriptors( I_j, labeled_elements_j );

d_i = [d_i(1:24,:); d_i2];
d_j = [d_j(1:24,:); d_j2];

n_bins_RGB_i = 24;
n_bins_RGB_j = 24;

n_bins_grad_i = 64;
n_bins_grad_j = 64;

SIFT_i = compute_SIFT_histograms( I_i, labeled_elements_i );
n_bins_SIFT_i = 128;
d_i = [d_i;SIFT_i];
SIFT_j = compute_SIFT_histograms( I_j, labeled_elements_j );
n_bins_SIFT_j = 128;
d_j = [d_j;SIFT_j];

%% Compute Q matrix (agirbau)
% Q matrix contains similarities between
% regions from both partitions
% fprintf(' Compute Q matrices...  \n')
%tic;
%Q11 = computeQIntra_sparse_mex(double(I_i),double(seg_i),idx_neighbors_i);
%Q22 = computeQIntra_sparse_mex(double(I_j),double(seg_j),idx_neighbors_j);

Q11_weights = compute_weights_QIntra_mex(idx_neighbors_i, ucm_i);
Q22_weights = compute_weights_QIntra_mex(idx_neighbors_j, ucm_j);

%Uncomment to combine default criteria (color and contour) with the UCM
%values criteria
% Q11 = Q11.*sparse(Q11_weights);
% Q22 = Q22.*sparse(Q22_weights);

%Uncomment to only use the UCM values criteria
Q11 = sparse(Q11_weights);
Q22 = sparse(Q22_weights);

[num_desc,~] = size(d_i);

Sigma = zeros(num_desc,num_desc);

for i=1:n_bins_RGB_i
    %Sigma(i,i) = 15;
    Sigma(i,i) = w_color;
end
for i=(n_bins_RGB_i+1):(n_bins_RGB_i+n_bins_grad_i)
    %Sigma(i,i) = 1;
    Sigma(i,i) = w_hog;
end
for i=(n_bins_RGB_i+n_bins_grad_i+1):(n_bins_RGB_i+n_bins_grad_i+n_bins_SIFT_i)
    Sigma(i,i) = w_sift;
end

W = computeWmatrix_sparse_motion_estimation_optical_flow_mex(d_i,d_j,labeled_elements_i,labeled_elements_j,Sigma,window_similarities, MV_x, MV_y);

%msg_trace = 'W matrix obtained'

Q12 = B_i'*W*B_j;   
    
Q= [Q11-1 Q12-1; Q12'-1 Q22-1];

% Q = [lambda*Q11 Q12; Q12' lambda*Q22]; 
%  Q= [lambda*Q11-0.2 Q12; Q12' lambda*Q22-0.2]; 
% Q= [lambda*Q11-0.2 0*Q12; 0*Q12' lambda*Q22-0.2];
% Q = [lambda*(Q11-0.3) Q12; Q12' lambda*(Q22-0.3)];
% Q = [lambda*(Q11-0.1) Q12; Q12' lambda*(Q22-0.1)];
%toc;
% fprintf(' ...ended! \n')
%%%%%%%%%%%%%%%%%%%%

%% Compute constraints 
% Compute constraints used in the optimization
% fprintf(' Compute constraints...  \n')
%tic;

[ A, b ] = inequalities_sparse(regions_i, regions_j, variables_opt, regions_adjacency, regions_variables );
A = A';
% A=[];
% b=[];



if hierarchical_constraint_opt

    [ nodes_leaves_i, variables_intra_i ] = compute_nodes_leaves( seg_i, merging_sequence_i );
    [ nodes_leaves_j, variables_intra_j ] = compute_nodes_leaves( seg_j, merging_sequence_j );
    
    [total_equalities2_i, total_inequalities2_i, b_equalities2_i, b_inequalities2_i] = compute_constraints_intra_hierarchy_mex(double(seg_i), double(merging_sequence_i), double(nodes_leaves_i), double(RAG_i), double(regions_variables), double(variables_opt), double(0));
    [total_equalities2_j, total_inequalities2_j, b_equalities2_j, b_inequalities2_j] = compute_constraints_intra_hierarchy_mex(double(seg_j), double(merging_sequence_j), double(nodes_leaves_j), double(RAG_j), double(regions_variables), double(variables_opt), double(max(max(seg_i))));

    A = [A; sparse(total_inequalities2_i); sparse(total_inequalities2_j)];
    b = [b b_inequalities2_i b_inequalities2_j];

    Aeq = [sparse(total_equalities2_i); sparse(total_equalities2_j)];
    beq = [b_equalities2_i b_equalities2_j];

else

    Aeq=[];
    beq=[];
    
end

edge_ids = select_edges_from_hierarchy( merging_sequence_i, RAG_i, regions_i, merging_sequence_j, RAG_j, regions_j, regions_adjacency);
Aeq2 = zeros(1,variables_opt);
Aeq2(edge_ids) = 1;
Aeq = [Aeq;Aeq2];
beq = [beq num_clusters-2];

% [edge_ids_i, edge_ids_j] = select_edges_from_hierarchy_separated( merging_sequence_i, RAG_i, regions_i, merging_sequence_j, RAG_j, regions_j, regions_adjacency);
% Aeq2 = zeros(1,variables_opt);
% Aeq2(edge_ids_i) = 1;
% Aeq3 = zeros(1,variables_opt);
% Aeq3(edge_ids_j) = 1;
% Aeq = [Aeq;Aeq2;Aeq3];
% beq = [beq num_clusters-1 num_clusters-1];


% Aeq = zeros(1, variables_opt);
% beq = 0;
%toc;
% fprintf(' ...ended! \n')
%%%%%%%%%%%%%%%%%%%%

%% Create cost function 
% fprintf(' Create cost function...  \n')
%tic;
f = zeros(1, variables_opt);

for i=1:length(regions_adjacency);
    region_i = regions_adjacency(1,i);
    region_j = regions_adjacency(2,i);
    
        f(i) = Q(region_i,region_j)+Q(region_j,region_i);
end
%toc;
% fprintf(' ...ended! \n')
%%%%%%%%%%%%%%%%%%%%

%% Optimization process 
% fprintf(' Optimization process...  \n')
%tic;
x = cplexbilp(f',A,b',Aeq,beq');
%toc;
% fprintf(' ...ended! \n')
%%%%%%%%%%%%%%%%%%%%

%% Create clustered partitions
% Create the resulting partitions using 
% the co-clustering result
% fprintf(' Create clustered partitions...  \n')
%tic;
if numel(x)>0
    [partition_i_inter, partition_j_inter, image_i_inter, image_j_inter] = create_clustered_partitions_original2(seg_i, seg_j, regions_adjacency, x);
    [partition_i_intra, partition_j_intra] = create_clustered_partitions_original2_intra(seg_i, seg_j, regions_adjacency, x);
else
    partition_i_inter = [];
    partition_j_inter = [];
    image_i_inter = [];
    image_j_inter = [];
    partition_i_intra = [];
    partition_j_intra = [];
end
%toc;
% fprintf(' ...ended! \n')
%%%%%%%%%%%%%%%%%%%%

%% Test the result
% [ result_i ] = check_partition( seg_i, nodes_leaves_i, partition_i_intra );
% [ result_j ] = check_partition( seg_j, nodes_leaves_j, partition_j_intra );
%%%%%%%%%%%%%%%%%%%%
% 
% figure;imshow(overlay_contour(uint8(image_i_inter),partition_i_inter,[0 0 0]))
% figure;imshow(overlay_contour(uint8(image_j_inter),partition_j_inter,[0 0 0]))

% figure;imshow(overlay_contour(label2rgb(partition_i_intra), partition_i_intra, [0 0 0]));
% figure;imshow(overlay_contour(label2rgb(partition_j_intra), partition_j_intra, [0 0 0]));

end
