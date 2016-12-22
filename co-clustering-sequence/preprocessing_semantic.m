function preprocessing_semantic(sequence_name,resolution)

    addpath('aux');
    score_th = 15;
    propagate_best_rnn_scores_intra_scale(sequence_name, score_th, resolution)

end