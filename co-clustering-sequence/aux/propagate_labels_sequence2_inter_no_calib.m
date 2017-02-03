function [ partitions_vector_out, images_vector_out, last_label ] = propagate_labels_sequence2_inter_no_calib( coarse_partition_i, coarse_partition_j, I_i, I_j, I_new, fine_partition_j, fine_partition_new, ucm_j, ucm_new, merging_sequence_j, merging_sequence_new, vectors_i, vectors_j, vectors_new, lambda, window_similarities, w_color, w_hog, w_sift, pixel_based_MV, hierarchical_constraint_MV, hierarchical_constraint_opt, motion_estimation, last_label )


%% Create variables
N_frames = 3;
[s_x, s_y] = size(coarse_partition_i);
labeled_elements_vector = zeros(2*s_x + 1, 2*s_y + 1, N_frames); %Vector containing all the labeled_elements
region_partitions       = zeros(1, N_frames);

contour_elems_partition = zeros(1, N_frames);                    %Number of contour elements at each partition

cumulate_regions           = 0;                                  %Total number of regions
variables_total            = 0;                                  %Total number of variables
Q_vector_total             = [];
regions_adjacency_total    = [];                                 %Vector containing all the adjacencies of the sequence
partitions_adjacency_total = [];                                 %Vector that relates each index of regions_adjacency_total
                                                                 %to the partition from which it belongs

d_vector = [];
cumulate_variables = 0;

B_struct = struct('frame', [], 'B_str',[], 'd_str', [], 'merging_sequence', [], ...
    'RAG', [], 'regions_variables', [], 'vectors', [], 'idx_neighbors', [], 'image', []);
B_vector = B_struct;
for ii=2:N_frames
    B_vector(ii) = B_struct;
end

B_vector(1).image =  I_i;
B_vector(2).image =  I_j;
B_vector(3).image =  I_new;

%% Relabel partitions
partitions_vector       = zeros(s_x, s_y, N_frames);             
[ relabeled_partition, labels_new1 ] = relabel_partition( coarse_partition_i );
partitions_vector(:,:,1) = relabeled_partition;


%% Similarities INTRA previous frames
%% Frame t-2
partition = partitions_vector(:,:,1);
B_vector(1).vectors = vectors_i;
region_partitions(1) = max(max(partition));

%% Obtain labeled elements
[ labeled_elements, gridbmap, idx_neighbors, angles ] = obtain_contour_elements( partition, vectors_i );
labeled_elements_vector(:,:,1)  = labeled_elements;
B_vector(1).idx_neighbors = idx_neighbors;

%% Compute B matrices 
[ B, gPb_angles ] = computeBmatrix_sparse_mex( double(partition), double(angles), double(labeled_elements));
B_vector(1).B_str = B;

%% Compute descriptors
[d, n_bins_RGB, n_bins_grad] = calculate_descriptors_mex(double(I_i),gPb_angles,labeled_elements, double(idx_neighbors.matrix_min), double(idx_neighbors.matrix_max));
[ d_2 ] = compute_HOG_descriptors( I_i, labeled_elements);
d = [d(1:24,:); d_2];
n_bins_RGB = 24;
n_bins_grad = 64;
SIFT = compute_SIFT_histograms( I_i, labeled_elements );
n_bins_SIFT = 128;
d = [d;SIFT];

B_vector(1).d_str = d;
contour_elems_partition(1) = size(d,2);

%% Compute Q intra
Qii = computeQIntra_sparse_mex(double(I_i),double(partition),idx_neighbors);
Qii = Qii-1;
%% Compute RAG intra
RAG_intra = RAG_Intra_mex(idx_neighbors, double(partition));
B_vector(1).RAG = RAG_intra;

%% Compute Adjacent regions
variables_opt = sum(sum(RAG_intra))/2;
[regions_adjacency, regions_variables] = variables_reduction(double(RAG_intra), double(variables_opt));

B_vector(1).regions_variables = regions_variables + cumulate_variables;
cumulate_variables = cumulate_variables + variables_opt;

%% Compute similarities intra vector
Qii_vector = zeros(1,size(regions_adjacency,2));
for ii=1:size(regions_adjacency,2)
    region_i = regions_adjacency(1,ii);
    region_j = regions_adjacency(2,ii);
    Qii_vector(ii) = Qii(region_i, region_j) + Qii( region_j, region_i);
end

%% Update regions adjacency
variables_total = variables_total + size(regions_adjacency,2);

Q_vector_tmp                  = Q_vector_total;
Q_vector_total                = zeros(1, variables_total);

regions_adjacency_total_tmp     = regions_adjacency_total;
regions_adjacency_total         = zeros(2, variables_total);

partitions_adjacency_total_tmp  = partitions_adjacency_total;
partitions_adjacency_total      = zeros(2, variables_total);

for ii=1:variables_total                
    if(ii <= size(regions_adjacency_total_tmp,2))    

        Q_vector_total(ii)             = Q_vector_tmp(ii);

        regions_adjacency_total(1,ii)    = regions_adjacency_total_tmp(1,ii);
        regions_adjacency_total(2,ii)    = regions_adjacency_total_tmp(2,ii);

        partitions_adjacency_total(1,ii) = partitions_adjacency_total_tmp(1,ii);
        partitions_adjacency_total(2,ii) = partitions_adjacency_total_tmp(2,ii);

    else

        Q_vector_total(ii)             = Qii_vector(ii-size(regions_adjacency_total_tmp,2));

        regions_adjacency_total(1,ii) = regions_adjacency(1,ii-size(regions_adjacency_total_tmp,2));
        regions_adjacency_total(2,ii) = regions_adjacency(2,ii-size(regions_adjacency_total_tmp,2));

        partitions_adjacency_total(1,ii) = 1;
        partitions_adjacency_total(2,ii) = 1;
    end

end

%% Update regions adjacency

%% Frame t-1
B_vector(2).merging_sequence = merging_sequence_j;
max_label = max(max(fine_partition_j));
partition = fine_partition_j;
labels_new2 = zeros(1, max_label);

for i=1:s_x
    for j=1:s_y        
        %labels_new2(partition(i,j)) = coarse_partition_i(i,j); %Check if it should be coarse_partition_j
        %If it is ok (coarse_partition_i), should I consider the motion
        %estimation?
        labels_new2(partition(i,j)) = coarse_partition_j(i,j);
    end
end

partitions_vector(:,:,2) = partition;

B_vector(2).vectors = vectors_j;
region_partitions(2) = max(max(partition));

%% Obtain labeled elements
[ labeled_elements, gridbmap, idx_neighbors, angles ] = obtain_contour_elements( partition, vectors_j );
labeled_elements_vector(:,:,2)  = labeled_elements;
B_vector(2).idx_neighbors = idx_neighbors;

%% Compute B matrices 
[ B, gPb_angles ] = computeBmatrix_sparse_mex( double(partition), double(angles), double(labeled_elements));
B_vector(2).B_str = B;

%% Compute descriptors
[d, n_bins_RGB, n_bins_grad] = calculate_descriptors_mex(double(I_j),gPb_angles,labeled_elements,double(idx_neighbors.matrix_min), double(idx_neighbors.matrix_max));
[ d_2 ] = compute_HOG_descriptors( I_j, labeled_elements );
d = [d(1:24,:); d_2];
n_bins_RGB = 24;
n_bins_grad = 64;
SIFT = compute_SIFT_histograms( I_j, labeled_elements );
n_bins_SIFT = 128;
d = [d;SIFT];

B_vector(2).d_str = d;
contour_elems_partition(2) = size(d,2);

%% Compute Q intra
Qii = computeQIntra_sparse_mex(double(I_j),double(partition),idx_neighbors);
Qii = lambda*Qii-1;
%% Compute RAG intra
RAG_intra = RAG_Intra_mex(idx_neighbors, double(partition));
B_vector(2).RAG = RAG_intra;

%% Compute Adjacent regions
variables_opt = sum(sum(RAG_intra))/2;
[regions_adjacency, regions_variables] = variables_reduction(double(RAG_intra), double(variables_opt));

B_vector(2).regions_variables = regions_variables + cumulate_variables;
cumulate_variables = cumulate_variables + variables_opt;

%% Compute similarities intra vector
Qii_vector = zeros(1,size(regions_adjacency,2));
for ii=1:size(regions_adjacency,2)
    region_i = regions_adjacency(1,ii);
    region_j = regions_adjacency(2,ii);
    Qii_vector(ii) = Qii(region_i, region_j) + Qii( region_j, region_i);
end

%% Update regions adjacency
variables_total = variables_total + size(regions_adjacency,2);

Q_vector_tmp                  = Q_vector_total;
Q_vector_total                = zeros(1, variables_total);

regions_adjacency_total_tmp     = regions_adjacency_total;
regions_adjacency_total         = zeros(2, variables_total);

partitions_adjacency_total_tmp  = partitions_adjacency_total;
partitions_adjacency_total      = zeros(2, variables_total);

for ii=1:variables_total                
    if(ii <= size(regions_adjacency_total_tmp,2))    

        Q_vector_total(ii)             = Q_vector_tmp(ii);

        regions_adjacency_total(1,ii)    = regions_adjacency_total_tmp(1,ii);
        regions_adjacency_total(2,ii)    = regions_adjacency_total_tmp(2,ii);

        partitions_adjacency_total(1,ii) = partitions_adjacency_total_tmp(1,ii);
        partitions_adjacency_total(2,ii) = partitions_adjacency_total_tmp(2,ii);

    else

        Q_vector_total(ii)             = Qii_vector(ii-size(regions_adjacency_total_tmp,2));

        regions_adjacency_total(1,ii) = regions_adjacency(1,ii-size(regions_adjacency_total_tmp,2));
        regions_adjacency_total(2,ii) = regions_adjacency(2,ii-size(regions_adjacency_total_tmp,2));

        partitions_adjacency_total(1,ii) = 2;
        partitions_adjacency_total(2,ii) = 2;
    end

end

%% Update regions adjacency


%% Similarities INTRA current frame
B_vector(3).merging_sequence = merging_sequence_new;
partition = fine_partition_new;
region_partitions(3)            = max(max(partition));
partitions_vector(:,:,3)        = partition;

%% Obtain labeled elements
[ labeled_elements, gridbmap, idx_neighbors, angles ] = obtain_contour_elements( partition, vectors_new );
labeled_elements_vector(:,:,3)  = labeled_elements;
B_vector(3).idx_neighbors = idx_neighbors;

%% Compute B matrices 
[ B, gPb_angles ] = computeBmatrix_sparse_mex( double(partition), angles, labeled_elements);
B_vector(3).B_str = B;

%% Compute descriptors
[d, n_bins_RGB, n_bins_grad] = calculate_descriptors_mex(double(I_new),gPb_angles,labeled_elements,double(idx_neighbors.matrix_min), double(idx_neighbors.matrix_max));
[ d_2 ] = compute_HOG_descriptors( I_new, labeled_elements );
d = [d(1:24,:); d_2];
n_bins_RGB = 24;
n_bins_grad = 64;
SIFT = compute_SIFT_histograms( I_new, labeled_elements );
n_bins_SIFT = 128;
d = [d;SIFT];

B_vector(3).d_str = d;
contour_elems_partition(3) = size(d,2);

%% Compute Q intra
Qii = computeQIntra_sparse_mex(double(I_new),double(partition),idx_neighbors);
Qii = lambda*Qii-1;
%% Compute RAG intra
RAG_intra = RAG_Intra_mex(idx_neighbors, double(partition));
B_vector(3).RAG = RAG_intra;

%% Compute Adjacent regions
variables_opt = sum(sum(RAG_intra))/2;
[regions_adjacency, regions_variables] = variables_reduction(double(RAG_intra), double(variables_opt));

 B_vector(3).regions_variables = regions_variables + cumulate_variables;
 cumulate_variables = cumulate_variables + variables_opt;

%% Compute similarities intra vector
Qii_vector = zeros(1,size(regions_adjacency,2));
for ii=1:size(regions_adjacency,2)
    region_i = regions_adjacency(1,ii);
    region_j = regions_adjacency(2,ii);
    Qii_vector(ii) = Qii(region_i, region_j) + Qii( region_j, region_i);
end   

%% Update regions adjacency
variables_total = variables_total + size(regions_adjacency,2);

Q_vector_tmp                  = Q_vector_total;
Q_vector_total                = zeros(1, variables_total);

regions_adjacency_total_tmp     = regions_adjacency_total;
regions_adjacency_total         = zeros(2, variables_total);

partitions_adjacency_total_tmp  = partitions_adjacency_total;
partitions_adjacency_total      = zeros(2, variables_total);

for ii=1:variables_total                
    if(ii <= size(regions_adjacency_total_tmp,2))    

        Q_vector_total(ii)             = Q_vector_tmp(ii);

        regions_adjacency_total(1,ii)    = regions_adjacency_total_tmp(1,ii);
        regions_adjacency_total(2,ii)    = regions_adjacency_total_tmp(2,ii);

        partitions_adjacency_total(1,ii) = partitions_adjacency_total_tmp(1,ii);
        partitions_adjacency_total(2,ii) = partitions_adjacency_total_tmp(2,ii);

    else

        Q_vector_total(ii)             = Qii_vector(ii-size(regions_adjacency_total_tmp,2));

        regions_adjacency_total(1,ii) = regions_adjacency(1,ii-size(regions_adjacency_total_tmp,2));
        regions_adjacency_total(2,ii) = regions_adjacency(2,ii-size(regions_adjacency_total_tmp,2));

        partitions_adjacency_total(1,ii) = 3;
        partitions_adjacency_total(2,ii) = 3;
    end

end
   
% variables_intra = length(Q_vector_total);

clearvars partition;
clearvars labeled_elements gridbmap idx_neighbors angles;
clearvars B gPb_angles;
clearvars d;
clearvars Q11 RAG_intra regions_adjacency regions_variables;
clearvars regions_adjacency_total_tmp partitions_adjacency_total_tmp;


%% Similarities INTER 
for ff=1:N_frames-1
    inter_relations = [ ff ff+1 ];
    if(ff > 1)
        inter_relations =[inter_relations;
                          ff-1 ff+1];
    end   
    
    for ll=1:size(inter_relations,1)
        ii = inter_relations(ll,1);
        jj = inter_relations(ll,2);

        fprintf(' Similarities inter: %d, %d  \n', ii,jj);
        %% Get descriptors 
        d_i = B_vector(ii).d_str;
        d_j = B_vector(jj).d_str;

        [num_desc,~] = size(d_i);

        %% Set variance
        Sigma = zeros(num_desc,num_desc);

        for i=1:n_bins_RGB
            Sigma(i,i) = w_color;
        end
        for i=(n_bins_RGB+1):(n_bins_RGB+n_bins_grad)
            Sigma(i,i) = w_hog;
        end
        for i=(n_bins_RGB+n_bins_grad+1):(n_bins_RGB+n_bins_grad+n_bins_SIFT)
            Sigma(i,i) = w_sift;
        end
        
 
        if motion_estimation
            [f_sift_i, d_sift_i] = vl_sift(im2single(rgb2gray(B_vector(ii).image)));
            [f_sift_j, d_sift_j] = vl_sift(im2single(rgb2gray(B_vector(jj).image)));
            [matches, scores] = vl_ubcmatch(d_sift_j, d_sift_i);
        %     matches = filter_matches(f_sift_j,f_sift_i,matches);
            matches = my_filter_matches(f_sift_j,f_sift_i,matches);

            if pixel_based_MV
                if hierarchical_constraint_MV
                    [idx_mat,estimated_MVs] = sift_based_motion_estimation_mat_hierarchical(partitions_vector(:,:,jj),B_vector(jj).merging_sequence,f_sift_j,f_sift_i,matches);
                else
                    [idx_mat,estimated_MVs] = sift_based_motion_estimation_mat(partitions_vector(:,:,jj),f_sift_j,f_sift_i,matches);
                end
                RAG_inter_ii_jj = RAG_Inter_motion_estimation_mex_pixel_based(double(partitions_vector(:,:,ii)), double(partitions_vector(:,:,jj)), idx_mat, estimated_MVs);
                W = computeWmatrix_sparse_motion_estimation_pixel_based_mex(d_i,d_j,double(labeled_elements_vector(:,:,ii)),double(labeled_elements_vector(:,:,jj)),double(Sigma),double(window_similarities), idx_mat, estimated_MVs);
            else
                if hierarchical_constraint_MV
                    motion_estimation_vectors = sift_based_motion_estimation_1vector_hierarchical(partitions_vector(:,:,jj),region_partitions(jj),B_vector(jj).merging_sequence,f_sift_j,f_sift_i,matches);
                else
                    motion_estimation_vectors = sift_based_motion_estimation_1vector(partitions_vector(:,:,jj),region_partitions(jj),f_sift_j,f_sift_i,matches);
                end
                RAG_inter_ii_jj = RAG_Inter_motion_estimation_mex(double(partitions_vector(:,:,ii)), double(partitions_vector(:,:,jj)), motion_estimation_vectors);
                W = computeWmatrix_sparse_motion_estimation_mex(d_i,d_j,double(labeled_elements_vector(:,:,ii)),double(labeled_elements_vector(:,:,jj)),double(Sigma),double(window_similarities), B_vector(jj).idx_neighbors, motion_estimation_vectors);
            end
        else
            RAG_inter_ii_jj = RAG_Inter_mex(double(partitions_vector(:,:,ii)), double(partitions_vector(:,:,jj)));
            W = computeWmatrix_sparse_mex(d_i,d_j,double(labeled_elements_vector(:,:,ii)),double(labeled_elements_vector(:,:,jj)),double(Sigma),double(window_similarities));
        end
        
        Qij = (B_vector(ii).B_str)'*W*(B_vector(jj).B_str);   


        %% Compute Adjacent regions
        variables_opt = sum(sum(RAG_inter_ii_jj));
        [regions_adjacency, regions_variables] = variables_reduction_inter(double(RAG_inter_ii_jj), double(variables_opt));

        %% Compute similarities intra vector
        if(ii==1)
            lambda_intra=3;
        else
            lambda_intra=1;
        end
        Qij_vector = zeros(1,length(regions_adjacency));
        for kk=1:size(regions_adjacency,2)
            Qij_vector(kk) = lambda_intra*(Qij(regions_adjacency(1,kk), regions_adjacency(2,kk))-1) ;
        end 

        %% Update regions adjacency
        variables_total = variables_total + size(regions_adjacency,2);

        Q_vector_tmp                  = Q_vector_total;
        Q_vector_total                = zeros(1, variables_total);

        regions_adjacency_total_tmp     = regions_adjacency_total;
        regions_adjacency_total         = zeros(2, variables_total);

        partitions_adjacency_total_tmp  = partitions_adjacency_total;
        partitions_adjacency_total      = zeros(2, variables_total);

        for kk=1:variables_total                
            if(kk <= length(regions_adjacency_total_tmp))    

                Q_vector_total(kk)             = Q_vector_tmp(kk);

                regions_adjacency_total(1,kk) = regions_adjacency_total_tmp(1,kk);
                regions_adjacency_total(2,kk) = regions_adjacency_total_tmp(2,kk);

                partitions_adjacency_total(1,kk) = partitions_adjacency_total_tmp(1,kk);
                partitions_adjacency_total(2,kk) = partitions_adjacency_total_tmp(2,kk);

            else

                Q_vector_total(kk)             = 2*real(Qij_vector(kk-length(regions_adjacency_total_tmp)));

                regions_adjacency_total(1,kk) = regions_adjacency(1,kk-length(regions_adjacency_total_tmp));
                regions_adjacency_total(2,kk) = regions_adjacency(2,kk-length(regions_adjacency_total_tmp));

                partitions_adjacency_total(1,kk) = ii;
                partitions_adjacency_total(2,kk) = jj;
            end

        end
    end
end

%% Conditions
[regions_adjacency_total_offset, regions_variables_offset, cumulate] = prepare_optimization(regions_adjacency_total, partitions_adjacency_total, region_partitions);

%% Compute constraints 
% Compute constraints used in the optimization
% fprintf(' Compute constraints...  \n')
%tic;
[ A, b ] = inequalities_sparse3(double(cumulate(end)), double(variables_total), double(regions_adjacency_total_offset), double(regions_variables_offset), double(partitions_adjacency_total) );
% [ A, b ] = inequalities_sparse(region_partitions(1), region_partitions(2), variables_total, regions_adjacency_total_offset, regions_variables_offset );
A = A';

Aeq = [];
beq = [];

%% Conditions INTRA 3
% for ii=3:3
%     [ nodes_leaves, variables_intra ] = compute_nodes_leaves( partitions_vector(:,:,ii), B_vector(ii).merging_sequence );
%     [total_equalities, total_inequalities, b_equalities, b_inequalities] = compute_constraints_intra_hierarchy_mex(double(partitions_vector(:,:,ii)), double(B_vector(ii).merging_sequence), double(nodes_leaves), double(B_vector(ii).RAG), double(regions_variables_offset), double(variables_total), double(cumulate(ii)));
% 
%     A = [A; sparse(total_inequalities)];
%     b = [b b_inequalities];
% 
%     Aeq = [Aeq; sparse(total_equalities)];
%     beq = [beq b_equalities];
% end

%% Conditions INTRA 1
%%Hay que mapear la informacion del vector x!!!!! 
condition_intra1 = zeros(1, length(partitions_adjacency_total));
count = 0;
for ii=1:length(partitions_adjacency_total)
    partition1 = partitions_adjacency_total(1,ii);
    partition2 = partitions_adjacency_total(2,ii);
    
    if(partition1==1 && partition2==1)
        condition_intra1(ii) = 1;
        count = count + 1;
    end
end
Aeq = [Aeq; condition_intra1];
beq = [beq count];

%% Conditions INTRA 2
%%Hay que mapear la informacion del vector x!!!!! 
condition_intra2_unions = zeros(1, length(partitions_adjacency_total));
condition_intra2_separates = zeros(1, length(partitions_adjacency_total));
count_separates = 0;
for ii=1:length(partitions_adjacency_total)
    partition1 = partitions_adjacency_total(1,ii);
    partition2 = partitions_adjacency_total(2,ii);
    
    if(partition1==2 && partition2==2)
        region1 = regions_adjacency_total(1,ii);
        region2 = regions_adjacency_total(2,ii);
        
        p = coarse_partition_j;
        
        index = p(partitions_vector(:,:,2)==region1);
        label1 = sum(index)/length(index);
        
        index = p(partitions_vector(:,:,2)==region2);
        label2 = sum(index)/length(index);
        
        if(label1==label2)
            condition_intra2_unions(ii) = 1;
        else
            condition_intra2_separates(ii) = 1;
            count_separates = count_separates +1;
        end
        
    end
end
Aeq = [Aeq; condition_intra2_unions];
beq = [beq 0];

Aeq = [Aeq; condition_intra2_separates];
beq = [beq count_separates];

%% Conditions INTER 12
unions = zeros(1, length(partitions_adjacency_total));
separates = zeros(1, length(partitions_adjacency_total));
count_separates = 0;

for ii=1:length(partitions_adjacency_total)
    partition1 = partitions_adjacency_total(1,ii);
    partition2 = partitions_adjacency_total(2,ii);
    
    if((partition1==1 && partition2==2)||(partition1==2 && partition1==1))        
        if(partition1==1)
            region1 = regions_adjacency_total(1,ii);
            region2 = regions_adjacency_total(2,ii);
        else
            region1 = regions_adjacency_total(2,ii);
            region2 = regions_adjacency_total(1,ii);
        end
        
        label1 = labels_new1(region1);
        label2 = region2;
        
        p2 = coarse_partition_j;
        index = p2(partitions_vector(:,:,2)==label2);
        label2 = sum(index)/length(index);
 
        if(label1==label2)
            unions(ii) = 1;
        else
            separates(ii) = 1;
            count_separates = count_separates +1;
        end
    end
end

Aeq = [Aeq; unions];
beq = [beq 0];

Aeq = [Aeq; separates];
beq = [beq count_separates];


tic;
addpath('/usr/local/opt/CPLEX_Studio124/cplex/matlab/')
x = cplexbilp(Q_vector_total',A,b',Aeq,beq');
toc;

if(~isempty(x))
% fprintf(' ...ended! \n')
%%%%%%%%%%%%%%%%%%%%
% [partition_i_inter, partition_j_inter, image_i_inter, image_j_inter] = create_clustered_partitions_original2(partitions_vector(:,:,1), partitions_vector(:,:,2), regions_adjacency_total_offset, x);
% 
% figure;imshow(overlay_contour(uint8(image_i_inter),partition_i_inter,[0 0 0]))
% figure;imshow(overlay_contour(uint8(image_j_inter),partition_j_inter,[0 0 0]))

% figure;imshow(overlay_contour(label2rgb(partition_i_intra), partition_i_intra, [0 0 0]));
% figure;imshow(overlay_contour(label2rgb(partition_j_intra), partition_j_intra, [0 0 0]));


partitions_vector_out(:,:,1) = coarse_partition_i;
partitions_vector_out(:,:,2) = coarse_partition_j;
% [partitions_vector_out] = create_clustered_partitions_sequence_intra(partitions_vector, regions_adjacency_total, partitions_adjacency_total, x);
[partitions_vector_out, images_vector_out, last_label] = create_clustered_partitions_sequence_propagate(partitions_vector, partitions_vector_out, regions_adjacency_total_offset, cumulate, x, labels_new1, labels_new2, last_label);
% [partitions_vector_out, images_vector_out] = create_clustered_partitions_sequence(partitions_vector, regions_adjacency_total_offset, cumulate, x);

% for i=1:N_frames
%     figure;imshow(overlay_contour(uint8(images_vector_out(:,:,:,i)),partitions_vector_out(:,:,i),[0 0 0]))
% %     path_out = sprintf('/work/dvaras/test2/partition_%03d.png', i);
% %     imwrite(overlay_contour(uint8(images_vector_out(:,:,:,i)),partitions_vector_out(:,:,i),[0 0 0]), path_out);
% end

else
partitions_vector_out = [];
images_vector_out = [];
last_label = [];
end
    

end