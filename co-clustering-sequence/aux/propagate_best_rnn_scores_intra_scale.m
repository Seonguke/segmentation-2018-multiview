function propagate_best_rnn_scores_intra_scale(sequence_name, score_th, resolution)

if ischar(score_th)
    score_th = str2num(score_th);
end

%Load partitions resulting from the co-clustering
load(['output/' sequence_name '/' sprintf('%03d',resolution) '_partition.mat']); 

%Load scores resulting from CRF as RNN
partitions_vector_total = partitions_vector_total_intra;
numframes = size(partitions_vector_total,3);
scores_vector = zeros(21,size(partitions_vector_total,1),size(partitions_vector_total,2),numframes);
sequence_path = 'data/images/';
frames = dir(strcat(fullfile(sequence_path,sequence_name),'/*.jpg'));

for ii=1:numframes
    image_filename_i = fullfile(sequence_path,sequence_name,frames(ii).name);
    [path_img_file, img_basename, img_extension] = fileparts(image_filename_i);
    load(['data/rnnscores/' sequence_name '/' img_basename '.mat']);
    for jj=1:21
        scores_label = reshape(scores(jj,:,:),size(scores,2),size(scores,3));
        scores_resized = imresize(scores_label,[size(partitions_vector_total,1) size(partitions_vector_total,2)]);
        scores_vector(jj,:,:,ii) = scores_resized;
    end
end

%Select the best scores and propagate them
classes_found_total = 0;
for ii=1:numframes
    partition = partitions_vector_total(:,:,ii);
    labels = unique(partition);
    for jj=1:numel(labels)
        region = (partition == labels(jj));
        stats_region = regionprops(region,'Area');
        area_region = 0;
        for mm=1:numel(stats_region)
            area_region = area_region + stats_region(mm).Area;
        end
        classes_found = 0;
        for kk=1:21
            scores_class = scores_vector(kk,:,:,ii);
            scores_class = reshape(scores_class,size(partition,1),size(partition,2));
            scores_class_region = scores_class(region);
            scores_class_region_above_th = (scores_class_region > score_th);
            stats_class = regionprops(scores_class_region_above_th,'Area');
            if ~isempty(stats_class)
                area_region_class = 0;
                for mm=1:numel(stats_class)
                    area_region_class = area_region_class + stats_class(mm).Area;
                end
                if area_region_class/area_region > 0.7
                    winning_classes(classes_found+1) = kk;
                    mean_values(classes_found+1) = mean(scores_class_region(scores_class_region_above_th));
                    classes_found = classes_found + 1;
                end
            end
        end
        if classes_found > 1
            sprintf('There are %02d winning classes in region %03d from frame %02d', classes_found, labels(jj), ii)
        end
        for kk=1:classes_found
            best_region_scores(classes_found_total+kk,:) = [mean_values(kk) winning_classes(kk) labels(jj) ii];
        end
        classes_found_total = classes_found_total + classes_found;
    end
end
[Y,I]=sort(best_region_scores(:,1),'descend');
best_regions_scores_sorted = best_region_scores(I,:);
output_path = ['output_semantic/' sequence_name]
mkdir(output_path);
save([output_path '/best_regions_intra_scale_' sprintf('%03d',resolution) '.mat'],'best_regions_scores_sorted');

end