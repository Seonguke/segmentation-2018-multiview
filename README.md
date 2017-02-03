# Semantically driven multiscale co-clustering for uncalibrated multiview segmentation

This code provides three different architectures for multiview uncalibrated segmentation:
- One-step iterative co-clustering (classical architecture)
- Two-step iterative co-clustering
- Two-step iterative co-clustering followed by a global optimization

The global optimization can take into account additional semantic information or be only based on the same general low-level features as those used in the iterative co-clustering.

Any of the previous architectures can be used either assuming that no motion is required or estimating the motion between the frames of the sequence.

## Pre-computing independent hierarchical segmentations

sequence_path: Base directory where sequence frames are stored

sequence_name: Name of the sequence

The frames of a sequence should be stored in sequence_path/sequence_name/ with JPG format.

ucm_dir: Base directory where UCMs (hierarchical segmentations) are stored

The UCMs of a sequence should be stored in ucm_dir/sequence_name/mat/ with MAT format. The basename of the UCM file, i.e. without the extension, should match the basename of the image filename. Furthermore, the probability contour values at each of the 8 orientations at which gPb is computed should be also stored in ucm_dir/sequence_name/vectors/ with MAT format and the filenames should also match the image basenames.

## Pre-computing independent semantic segmentations (only for semantic segmentation)

## Computing co-clustering
