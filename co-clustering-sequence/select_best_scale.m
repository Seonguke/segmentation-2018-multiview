function select_best_scale(sequence_name, output_path, results_path)

addpath('aux');

%Load partitions resulting from the co-clustering
partition_files = dir([results_path '/' sequence_name '/*_partition.mat']);

num_scales = numel(partition_files);
for ii=1:num_scales  
    load([results_path '/' sequence_name '/' partition_files(ii).name]);    
    partitions(:,:,:,ii) = partitions_vector_total;
end

%Load CNN scores
numframes = size(partitions,3);
scores_vector = zeros(21,size(partitions_vector_total,1),size(partitions_vector_total,2),numframes);
sequence_path = 'data/images/';
frames = dir(strcat(fullfile(sequence_path,sequence_name),'/*.jpg'));
for ii=1:numframes
    image_filename = fullfile(sequence_path,sequence_name,frames(ii).name);
    [path_img_file, img_basename, img_extension] = fileparts(image_filename);
    load(['data/rnnscores/' sequence_name '/' img_basename '.mat']);
    for jj=1:21
        scores_label = reshape(scores(jj,:,:),size(scores,2),size(scores,3));
        scores_resized = imresize(scores_label,[size(partitions_vector_total,1) size(partitions_vector_total,2)]);
        scores_vector(jj,:,:,ii) = scores_resized;
    end
end

result_file = [output_path '/' sequence_name '/clusters_selected.mat'];


if ~exist(result_file)

    %Compute semantic segmentation outputs from CNN scores
    for ii=1:numframes
        scores_ii = scores_vector(:,:,:,ii);
        [M,I] = max(scores_ii,[],1);
        semantic_segmentation(:,:,ii) = reshape(I,size(scores_ii,2),size(scores_ii,3));
    end

    J_scales = zeros(1,num_scales);
    clusters_selected = [];
    min_valid_label = 1;
    
    for ii=1:num_scales
        J_seq = 0;
        for jj=1:numframes
            J_frame = 0;
            labels = unique(semantic_segmentation(:,:,jj));
            num_labels = 0;
            for kk=1:numel(labels)
                label = labels(kk);
                if label > min_valid_label
                    num_labels = num_labels + 1;
                    semantic_mask = semantic_segmentation(:,:,jj)==label;
                    [J, part_selected] = partition_upper_bound_J(partitions(:,:,jj,ii),semantic_mask);
                    for mm=2:numel(unique(part_selected))
                        [row,col] = find(part_selected == mm, 1, 'first');
                        label_cluster = partitions(row,col,jj,ii);
                        clusters_selected = [clusters_selected; ii jj label label_cluster];
                    end
                    J_frame = J_frame + J;
                end
            end
            J_frame = J_frame / num_labels;
            J_seq = J_seq + J_frame;
        end
        J_scales(ii) = J_seq/numframes;
    end

    save(result_file,'clusters_selected');
else
    load(result_file);
end



for ii=1:21
    scores_appended_label = [];
    for jj=1:numframes
        scores_tmp = scores_vector(ii,:,:,jj);
        scores_tmp = reshape(scores_tmp,size(scores_tmp,2),size(scores_tmp,3));
        scores_appended_label = [scores_appended_label scores_tmp];
    end
    scores_appended{ii} = scores_appended_label;
end

image_filename = fullfile(sequence_path,sequence_name,frames(1).name);
[path_img_file, img_basename, img_extension] = fileparts(image_filename);
[im_drop,color_map] = imread(['data/semantic_segs/' sequence_name '/' img_basename '.png']);

for ii=1:num_scales
    %Compute scores for each cluster and each semantic label
    
    partitions_appended = [];
    for jj=1:numframes
        partitions_appended = [partitions_appended partitions(:,:,jj,ii)];
    end
    
    rows = find(clusters_selected(:,1)==ii);
    clusters_scale = clusters_selected(rows,:);
    
    scores_table = [];
    for jj=1:21
        scores_rnn_label = scores_appended{jj};
        rows = find(clusters_scale(:,3)==jj);
        clusters_semantic_label = clusters_scale(rows,4);
        clusters_semantic_label = unique(clusters_semantic_label);
        for kk=1:numel(clusters_semantic_label)
            mask = (partitions_appended==clusters_semantic_label(kk));
            scores_cluster = sum(scores_rnn_label(mask));
            if isempty(scores_table)
                scores_table = [clusters_semantic_label(kk) scores_cluster jj];
            else
                pos = find(scores_table(:,1)==clusters_semantic_label(kk));
                if isempty(pos)
                    scores_table = [scores_table; clusters_semantic_label(kk) scores_cluster jj];
                else
                    %Remove non-maximum semantic clusters when there are conflicts about
                    %semantic in the clusters
                    if scores_cluster > scores_table(pos,2)
                        scores_table(pos,2) = scores_cluster;
                        scores_table(pos,3) = jj;
                    end
                end
            end
        end
    end
            
    %Compute score for each scale and show the resulting semantic
    %segmentation for each scale
    score_scale{ii} = sum(scores_table(:,2));
    semantic_partition_appended = zeros(size(partitions_appended));
    mask_selected_clusters_appended = false(size(partitions_appended));
    for jj=1:21
        rows = find(scores_table(:,3)==jj);
        for kk=1:numel(rows)
            cluster = scores_table(rows(kk),1);
            semantic_partition_appended(partitions_appended==cluster)=jj;
            mask_selected_clusters_appended(partitions_appended==cluster)=true;
        end
    end
    semantic_partitions{ii} = semantic_partition_appended;
    %figure;imshow(semantic_partition_appended,color_map);
    
    mask_background_appended = ~mask_selected_clusters_appended;
    scores_rnn_bg = scores_appended{1};
    %score_scale{ii}
    score_scale{ii} = score_scale{ii} + sum(scores_rnn_bg(mask_background_appended));
    %[sum(scores_rnn_bg(mask_background_appended)) score_scale{ii}]
    
end

%save(['semantic_segmentations_' sequence_name '.mat'], 'semantic_partitions');

%Select the scale with best score
max_score = -1;
best_scale = -1;
for ii=1:num_scales
    if score_scale{ii} > max_score
        max_score = score_scale{ii};
        best_scale = ii;
    end
end
best_semantic_segmentation = semantic_partitions{best_scale};

output_file = [output_path '/' sequence_name '/best_scale.mat'];

save(output_file,'best_scale','max_score','best_semantic_segmentation','score_scale');

        
end