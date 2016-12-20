function co_clustering_sequence_multiview_no_calibration(sequence_name, out_path, num_regions)

if nargin == 2
    num_regions = 10;
end

sequence_path = 'data/images/';
ucm_dir = 'data/ucms/';


window = 20;
w_color = 15;
w_hog = 1;
w_sift = 0;
sparseSIFT = false;
init_reg = 200;
motion_estimation = false;
pixel_based_MV = false;
optical_flow = false;
hierarchical_constraint_MV = true;
hierarchical_constraint_opt = true;
set(0,'RecursionLimit',1500);

path_ima = sprintf('%s/%s/%03d_partition.mat', out_path, sequence_name, num_regions);
if ~exist(path_ima)
    
    %% Create output folders
    create_folder=sprintf('mkdir %s/%s', out_path, sequence_name);
    system(create_folder);

    give_acces=sprintf('chmod 777 %s/%s', out_path, sequence_name);
    system(give_acces);

    frames = dir(strcat(fullfile(sequence_path,sequence_name),'/undistorted*.jpg'));
    ucm_video_path = fullfile(ucm_dir,sequence_name,'mat');
    ucm_vectors_video_path = fullfile(ucm_dir,sequence_name,'vectors');


    image_filename_i = fullfile(sequence_path,sequence_name,frames(1).name);
    I_i = imread(image_filename_i);
    [path_img_file, img_basename, img_extension] = fileparts(image_filename_i);
    ucm_filename_i = fullfile(ucm_video_path,[img_basename '.mat']);
    ucm_i = load(ucm_filename_i);
    ucm_i = ucm_i.ucm2;
    vector_filename_i = fullfile(ucm_vectors_video_path,[img_basename '.mat']);
    vectors_i = load(vector_filename_i);
    vectors_i = vectors_i.ws_a;
    if init_reg~=0
        [ seg_i, merging_sequence_i ] = ucm2bpt(ucm_i,init_reg);
    else
        [ seg_i, merging_sequence_i ] = ucm2bpt(ucm_i);
    end

    image_filename_j = fullfile(sequence_path,sequence_name,frames(2).name);
    I_j = imread(image_filename_j);
    [path_img_file, img_basename, img_extension] = fileparts(image_filename_j);
    ucm_filename_j = fullfile(ucm_video_path,[img_basename '.mat']);
    ucm_j = load(ucm_filename_j);
    ucm_j = ucm_j.ucm2;
    vector_filename_j = fullfile(ucm_vectors_video_path,[img_basename '.mat']);
    vectors_j = load(vector_filename_j);
    vectors_j = vectors_j.ws_a;
    if init_reg~=0
        [ seg_j, merging_sequence_j ] = ucm2bpt(ucm_j,init_reg);
    else
        [ seg_j, merging_sequence_j ] = ucm2bpt(ucm_j);
    end

    [partition_i_inter, partition_j_inter, image_i_inter, image_j_inter, partition_i_intra, partition_j_intra] = ...
        co_clustering_global_hierarchy_num_clusters( I_i, I_j, seg_i, seg_j, ucm_i, ucm_j, merging_sequence_i, merging_sequence_j, vectors_i, vectors_j, window, w_color, w_hog, w_sift, pixel_based_MV, hierarchical_constraint_MV, hierarchical_constraint_opt, motion_estimation, 2*num_regions  );

    label1 = max(max(partition_i_inter));
    label2 = max(max(partition_j_inter));
    last_label=max([label1 label2]);

    partitions_vector_total(:,:,1) = partition_i_inter;
    partitions_vector_total(:,:,2) = partition_j_inter;

    for ii=3:numel(frames)
        ii
        image_filename_new = fullfile(sequence_path,sequence_name,frames(ii).name);
        I_new = imread(image_filename_new);
        [path_img_file, img_basename, img_extension] = fileparts(image_filename_new);
        ucm_filename_new = fullfile(ucm_video_path,[img_basename '.mat']);
        ucm_new = load(ucm_filename_new);
        ucm_new = ucm_new.ucm2;
        vector_filename_new = fullfile(ucm_vectors_video_path,[img_basename '.mat']);
        vectors_new = load(vector_filename_new);
        vectors_new = vectors_new.ws_a;
        if init_reg~=0
            [ seg_new, merging_sequence_new ] = ucm2bpt(ucm_new,init_reg);
        else
            [ seg_new, merging_sequence_new ] = ucm2bpt(ucm_new);
        end
        [ partitions_vector_out, images_vector_out, last_label ] = propagate_labels_sequence2_no_calibration( partition_i_inter, partition_j_inter, I_i, I_j, I_new, seg_j, seg_new, ucm_j, ucm_new, merging_sequence_j, merging_sequence_new, vectors_i, vectors_j, vectors_new, window, w_color, w_hog, w_sift, pixel_based_MV, hierarchical_constraint_MV, hierarchical_constraint_opt, motion_estimation, 2*num_regions, last_label);

        if ~isempty(partitions_vector_out)
            partitions_vector_total(:,:,ii) = partitions_vector_out(:,:,3);
            partition_i_inter = partitions_vector_out(:,:,2);
            partition_j_inter = partitions_vector_out(:,:,3);
            I_i = I_j;
            I_j = I_new;
            seg_j = seg_new;
            ucm_j = ucm_new;
            merging_sequence_j = merging_sequence_new;
            vectors_i = vectors_j;
            vectors_j = vectors_new;
        else
            break
        end
    end
    if ii == numel(frames)
        save(path_ima, 'partitions_vector_total');
    end
end

% max_label = max(max(max(partitions_vector_total)));
% tmp = jet(double(max_label));
% color = uint8(255*tmp(randperm(max_label,max_label),:));
% 
% images_vector_total = zeros(size(partitions_vector_total,1),size(partitions_vector_total,2),3,size(partitions_vector_total,3));
% for kk=1:size(images_vector_total,4)
%     for i=1:size(images_vector_total,1)
%         for j=1:size(images_vector_total,2)
%             images_vector_total(i,j,1,kk) = color(partitions_vector_total(i,j,kk),1);
%             images_vector_total(i,j,2,kk) = color(partitions_vector_total(i,j,kk),2);
%             images_vector_total(i,j,3,kk) = color(partitions_vector_total(i,j,kk),3);
%         end
%     end
%     figure;imshow(overlay_contour(uint8(images_vector_total(:,:,:,kk)),partitions_vector_total(:,:,kk),[0 0 0]));
% end

end
