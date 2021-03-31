# GISTIC2 postprocessing

To force-call CNV drivers for each tumor sample.


Requirements: Matlab

The inputs to the Matlab script are: segmentation_file, amp_focal_file (GISTIC2 output), del_focal_file (GISTIC2 output), cnv_blacklist_file (GISTIC2 input), arm_coordinates_file. 

The outputs contain 6 files:
1. presence/absence of amplification/deletion per chromosome arm per sample;
2. median intensity per chromosome arm per sample;
3. presence/absence of each focal amplification driver per sample; 
4. corrected intensity of each focal amplification driver per sample;
5. presence/absence of each focal deletion driver per sample;
6. corrected intensity of each focal deletion driver per sample;

Steps:

1. Open Matlab

2. Run the following commands to assign the filepaths to the required inputs:

clear

segfile = 'example_data/concat_seg.aggregated.exclude.nan.tsv';

amp_focal_file = 'example_data/amp_focal_file.txt';

del_focal_file = 'example_data/del_focal_file.txt';

cnv_blacklist_file = 'example_data/WES_pairs_rerun_final_20200703_update_blacklist.txt';

arm_coordinates_file = 'example_data/hg19.GISTIC.arms.tsv';

3. Then execute the function GISTIC2_force_calling:

GISTIC2_force_calling(segfile,amp_focal_file,del_focal_file,cnv_blacklist_file, arm_coordinates_file)