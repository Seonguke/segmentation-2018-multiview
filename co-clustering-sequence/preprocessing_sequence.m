function preprocessing_sequence( sequence_name )

addpath('aux');
addpath('../external/ucm/');
addpath('../external/optical_flow/');

video_dir = ['data/images/' sequence_name];
ucm_dir = ['data/ucms/' sequence_name];
mkdir(ucm_dir);

frames = dir(strcat(video_dir,'/*.jpg'));
ucm_video_path = fullfile(ucm_dir,'mat');
mkdir(ucm_video_path);
ucm_vectors_video_path = fullfile(ucm_dir,'vectors');
mkdir(ucm_vectors_video_path);

%% Computing UCM
for jj=1:numel(frames)
    image_filename = fullfile(video_dir,frames(jj).name);
    [path_img_file, img_basename, img_extension] = fileparts(image_filename);
    ucm_filename = fullfile(ucm_video_path,[img_basename '.mat']);
    [ucm, ws_a] = im2ucm_normal_vectors(image_filename,ucm_filename);
    vector_filename = fullfile(ucm_vectors_video_path,[img_basename '.mat']);
    save(vector_filename,'ws_a');
end

%% Computing optical flow between consecutive frames
optflow_dir = ['data/optical_flow/' sequence_name '/previous_frame/'];
mkdir(optflow_dir);
image_filename = fullfile(video_dir,frames(1).name);
previous_frame = imread(image_filename);
for jj=2:numel(frames)
    image_filename = fullfile(video_dir,frames(jj).name);
    current_frame = imread(image_filename);
    save_path = sprintf([optflow_dir '/frame_OF_%05d.mat'], jj);
    if ~exist(save_path)    
        flow = mex_LDOF(double(current_frame), double(previous_frame));
        save(save_path, 'flow');
    end    
    previous_frame = current_frame;
end
    
%% Computing optical flow between frame at distance 2
optflow_dir = ['data/optical_flow/' sequence_name '/2previous_frames/'];
mkdir(optflow_dir);
for jj=3:numel(frames)
    image_filename = fullfile(video_dir,frames(jj-2).name);
    previous_frame = imread(image_filename);
    image_filename = fullfile(video_dir,frames(jj).name);
    current_frame = imread(image_filename);
    save_path = sprintf([optflow_dir '/frame_OF_%05d.mat'], jj);
    if ~exist(save_path)    
        flow = mex_LDOF(double(current_frame), double(previous_frame));
        save(save_path, 'flow');
    end    
end

end

