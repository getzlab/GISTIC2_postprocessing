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

2. Run the following commands to assign the filepath to the required inputs:

clear

segfile = 'segmentation_file';

amp_focal_file = 'amp_focal_file';

del_focal_file = 'del_focal_file';

cnv_blacklist_file = 'cnv_blacklist_file';

arm_coordinates_file = 'arm_coordinates_file';

3. Then execute the function GISTIC2_force_calling:

GISTIC2_force_calling(segfile,amp_focal_file,del_focal_file,cnv_blacklist_file, arm_coordinates_file)