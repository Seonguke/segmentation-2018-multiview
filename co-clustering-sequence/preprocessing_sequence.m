function preprocessing_sequence( sequence_name )

addpath('aux');
addpath('../external/ucm/');

video_dir = ['data/images/' sequence_name];
ucm_dir = ['data/ucms/' sequence_name];
mkdir(ucm_dir);

frames = dir(strcat(video_dir,'/*.jpg'));
ucm_video_path = fullfile(ucm_dir,'mat');
mkdir(ucm_video_path);
ucm_vectors_video_path = fullfile(ucm_dir,'vectors');
mkdir(ucm_vectors_video_path);
for jj=1:numel(frames)
    image_filename = fullfile(video_dir,frames(jj).name);
    [path_img_file, img_basename, img_extension] = fileparts(image_filename);
    ucm_filename = fullfile(ucm_video_path,[img_basename '.mat']);
    [ucm, ws_a] = im2ucm_normal_vectors(image_filename,ucm_filename);
    vector_filename = fullfile(ucm_vectors_video_path,[img_basename '.mat']);
    save(vector_filename,'ws_a');
end


end

