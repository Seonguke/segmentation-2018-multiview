function co_clustering_global_after_iterative(sequence_name, out_path, results_path, num_regions, window_size)

if nargin == 3
    num_regions = 10;
    window_size = 2;
elseif nargin == 4
    window_size = 2;
end

addpath('aux');
addpath('aux/mex');
addpath(fullfile('../external/vlfeat-0.9.20','toolbox','sift')) ;
addpath(fullfile('../external/vlfeat-0.9.20','toolbox','mex','mexa64')) ;
addpath('../external/cplex/')

sequence_path = 'data/images/';
ucm_dir = 'data/ucms/';

window_similarities = 20;
w_color = 15;
w_hog = 1;
w_sift = 0;
set(0,'RecursionLimit',1500);

path_ima = sprintf('%s/%s/%03d_partition.mat', out_path, sequence_name, num_regions);
if ~exist(path_ima)
    
    mkdir([out_path '/' sequence_name]);
    
    path_inter = [results_path '/' sequence_name '/' sprintf('%03d',num_regions) '_partition.mat'];
    if exist(path_inter)
        prev_results = load(path_inter);
        partitions_vector_total_intra = prev_results.partitions_vector_total;
        
        frames = dir(strcat(fullfile(sequence_path,sequence_name),'/*.jpg'));
        ucm_video_path = fullfile(ucm_dir,sequence_name,'mat');
        ucm_vectors_video_path = fullfile(ucm_dir,sequence_name,'vectors');
            
        %Global optimization over intra partitions obtained by an iterative approach
            
        cumulate_regions           = 0;
        variables_total            = 0;
        Q_vector_total             = [];
        regions_adjacency_total    = [];
        partitions_adjacency_total = [];
        
        d_vector = [];
        cumulate_variables = 0;
        
        N_frames = numel(frames);
        if size(partitions_vector_total_intra,3) == N_frames   
            region_partitions       = zeros(1, N_frames);
            contour_elems_partition = zeros(1, N_frames); 
            image_filename = fullfile(sequence_path,sequence_name,frames(1).name);
            I = imread(image_filename);
            s_x = size(I,1);
            s_y = size(I,2);
            partitions_vector = zeros(s_x,s_y,N_frames);
            labeled_elements_vector = zeros(2*s_x + 1, 2*s_y + 1, N_frames);
            
            B_struct = struct('frame', [], 'B_str',[], 'd_str', [], 'merging_sequence', [], ...
            'RAG', [], 'regions_variables', [], 'vectors', [], 'idx_neighbors', [], 'image', []); 
            
            %% Similarities INTRA frames
            for ff=1:N_frames
                
                fprintf(' Similarities intra: %d\n', ff);
                B_vector(ff) = B_struct;
                
                image_filename = fullfile(sequence_path,sequence_name,frames(ff).name);
                I = imread(image_filename);
                B_vector(ff).image = I;
                [path_img_file, img_basename, img_extension] = fileparts(image_filename);

                vector_filename = fullfile(ucm_vectors_video_path,[img_basename '.mat']);
                vectors = load(vector_filename);
                vectors = vectors.ws_a;
                seg = relabel_partition(partitions_vector_total_intra(:,:,ff));
                merging_sequence = [];
            
                partitions_vector(:,:,ff) = seg;
                B_vector(ff).vectors = vectors;
                region_partitions(ff) = max(max(seg));
                B_vector(ff).merging_sequence = merging_sequence;
                
                %% Obtain labeled elements
                [ labeled_elements, gridbmap, idx_neighbors, angles ] = obtain_contour_elements( seg, vectors );
                labeled_elements_vector(:,:,ff)  = labeled_elements;
                B_vector(ff).idx_neighbors = idx_neighbors;
                
                %% Compute B matrices 
                [ B, gPb_angles ] = computeBmatrix_sparse_mex( double(seg), double(angles), double(labeled_elements));
                B_vector(ff).B_str = B;
                
                %% Compute descriptors
                [d, n_bins_RGB, n_bins_grad] = calculate_descriptors_mex(double(I),gPb_angles,labeled_elements, double(idx_neighbors.matrix_min), double(idx_neighbors.matrix_max));
                [ d_2 ] = compute_HOG_descriptors( I, labeled_elements);
                d = [d(1:24,:); d_2];
                n_bins_RGB = 24;
                n_bins_grad = 64;
                SIFT = compute_SIFT_histograms( I, labeled_elements );
                n_bins_SIFT = 128;
                d = [d;SIFT];
        
                B_vector(ff).d_str = d;
                contour_elems_partition(ff) = size(d,2);
                
                %% Compute Q intra
                Qii = computeQIntra_sparse_mex(double(I),double(seg),idx_neighbors);
                Qii = Qii-1;
                
                %% Compute RAG intra
                RAG_intra = RAG_Intra_mex(idx_neighbors, double(seg));
                B_vector(ff).RAG = RAG_intra;
        
                %% Compute Adjacent regions
                variables_opt = sum(sum(RAG_intra))/2;
                [regions_adjacency, regions_variables] = variables_reduction(double(RAG_intra), double(variables_opt));
                
                B_vector(ff).regions_variables = regions_variables + cumulate_variables;
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
        
                        partitions_adjacency_total(1,ii) = ff;
                        partitions_adjacency_total(2,ii) = ff;
                        
                    end
                end
            end        
            
            %% Similarities INTER frames
            for ff=1:N_frames-1
                inter_relations = [ ff ff+1 ];
                last_frame = min(N_frames,ff+window_size-1);
                for mm = ff+2:last_frame
                    inter_relations =[inter_relations;
                                    ff mm];
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
                    
                    if jj==ii+1
                        load(['data/optical_flow/' sequence_name '/previous_frame/frame_OF_' sprintf('%05d',jj) '.mat']);
                    else
                        diff_idxs = jj-ii;
                        load(['data/optical_flow/' sequence_name '/' sprintf('%d',diff_idxs) 'previous_frames/frame_OF_' sprintf('%05d',jj) '.mat']);
                    end
                    
                    MV_x = flow(:,:,1);
                    MV_y = flow(:,:,2);
                    RAG_inter_ii_jj = RAG_Inter_motion_estimation_optical_flow_mex(double(partitions_vector(:,:,ii)), double(partitions_vector(:,:,jj)), MV_x, MV_y);
                    W = computeWmatrix_sparse_motion_estimation_optical_flow_mex(d_i,d_j,double(labeled_elements_vector(:,:,ii)),double(labeled_elements_vector(:,:,jj)),double(Sigma),double(window_similarities), MV_x, MV_y);
        
                    Qij = (B_vector(ii).B_str)'*W*(B_vector(jj).B_str);
                    
                    %% Compute Adjacent regions
                    variables_opt = sum(sum(RAG_inter_ii_jj));
                    [regions_adjacency, regions_variables] = variables_reduction_inter(double(RAG_inter_ii_jj), double(variables_opt));
                    
                    %% Compute normalization vector
                    %normalization_vector = obtain_normalization_vector(regions_adjacency, B_vector(ii).idx_neighbors, B_vector(jj).idx_neighbors, region_partitions(ii), region_partitions(jj));
                    %max_val = max(normalization_vector);
                    
                    %% Compute similarities intra vector
                    Qij_vector = zeros(1,length(regions_adjacency));
                    for kk=1:size(regions_adjacency,2)
                        Qij_vector(kk) = Qij(regions_adjacency(1,kk), regions_adjacency(2,kk))-1;
                        %Qij_vector(kk) = 1/normalization_vector(kk)*Qij(regions_adjacency(1,kk), regions_adjacency(2,kk));
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
            [ A, b ] = inequalities_sparse3(double(cumulate(end)), double(variables_total), double(regions_adjacency_total_offset), double(regions_variables_offset), double(partitions_adjacency_total) );
            A = A';
            
            Aeq = [];
            beq = [];
                
            tic;
            x = cplexbilp(Q_vector_total',A,b',Aeq,beq');
            toc;   
            
            if(~isempty(x))
                [partitions_vector_total, images_vector_out] = create_clustered_partitions_sequence_global_opt(partitions_vector, regions_adjacency_total_offset, cumulate, x);
                save(path_ima, 'partitions_vector_total');
            end
        end
    end
end

end