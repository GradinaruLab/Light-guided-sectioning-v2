# Light-guided-sectioning-v2
Cell Reports 2021

Figure 1: 
Find features of slicing images
!. Take images at 50um steps, at 3 conditions- bright field, dark- low exp, dark- high exposure. Image the cleared brain to identify the distance from implant surface. 
Upload to computer. Separate to 3 folders. 
Run: AnalyzeStackOfImages – you will need SIFT and vlfeat-0.9.21 libraries to run it. Before you run 
That will create .mat files for each distance
D:\DATA_Glab\Light_sectioning\LiGS pictures cryosectioning\sample1\run_extract_features.m
Repeat that for all samples. 
this function look at the gaussian fits of samples that had implants with the same diameter:  ( first run AnalyseStackOfImages)

extract_gaussian_fit_features_LiGS_slicing_by_diameter

Figure 2: 
(and Figure S2) 
For FP SCN analysis: 
‘Test6Rcompare’- uses to look at ‘test6R’ exaperiments and also compare signal to LiGS histology, over samples. It reads information from ‘LGS_SCN_VIP’ xlsx file
Figure 3: 
Image quantification: 
C:\Users\anatk\Documents\Matlab\LightGuidedSectioning\get_fluor_quantification_over_samples_LGS_042820.m
Uses: 
find_fluor_quantification_zstack_10umStep_042820(mouse)
has also:  ‘count_cells’

Check the number of cells identified in in-vivo vs fixed (through GRIN):
‘run_2P_registration_cell_quantification’ which Uses: ‘LGS_2P_registration_cell_quantifiction’
also compare to shuffled cells. The basic function used is ‘ssim’
After registration, use: 
Run_extract_registration_parameters_LGS
To analyse the distances, angle and scale of the registration 

Figure 4: 

2P

1)	Image z-stack 2P tissue side, through GRIN lens, and confocal, tissue side, if more colors are available. Tissue side can be imaged from 0 (GRIN lens surface) to ~200/250 um, as more is not relevant. GRIN side- all the range that have signal also in the in-vivo imaging. 
2)	Invivo data analysis: 
In vivo, 3 min recording each Z plane : 
a.	to save tiff files for each 3 minutes of Z-stack I use ‘Zstack_show_movie_persection’ that show movies and save them (as movies). 
b.	Apply motion correction: apply_motion_correction(‘non-rigid’) or (‘rigid’)
c.	To save tiff files for suite 2P:
Y=read_file(‘file name’);
To take a look at a frame, for example f=1: imshow(unit16(Y(:,:,1))
Save using: saveastiff(uint16(Y),’new file name.tiff’)
d.	To apply motion correction and cell identification, suite2P (python version) is used, using the following parameters: tau=1.25, fs=4, save_mat=1, nonrigid=1, block size= 128,128, sparse_mode=1. The rest of parameters remains the default parameters.  
In cases where suite2P failed to produce a good motion correction, the function ‘apply_motion_correction’ is used, which apply motion correction to all the files in the folder, either rigid or non-rigid (rigid_MC / nonrigid_MC, written by Ryan Cho). Also save new MEAN images. Cells activity was inspected using suite2P afterwards.  

3)	2p-fixed tissue Zsatck : create a new folder for the new .dat. While this folder is open (empty), run ‘plot_all_Zstack_2P’ and choose the relevant .dat file 

4)	Inspect Z-stack using imageJ. Find few clear marks, that following them allow a full reconstruction of all z stack. 
5)	Using ‘Register_2D_Zstack_2P_v2’ – register selected planes for each brain. For each transformation saves correlation coefficient (R), transform matrix (tform), and angle and scale. This function also uses 'LGS_get_angle_and_scale' to get the angle and the scale of the transformation for further analysis. 	
v2 has the option to register directly from in-vivo to Tissue, as I needed for WT2170 and WT2174
7) Optional: Volume registration: ‘make_complete_zstack_registration’
8)	Apply_transformation_on_stained_image – apply transformation on 647 staining and save it

Figure 5 and 6: 

In vivo to stained images: 
8)	To register suite2P data to stained images- 
‘get_data_suite2P’- which creates the contour maps of all cells, of stained cells (if there are) and creates folder of images that can be uploads using imageJ for comparing in-vivo to stained image.*** 
‘run_get_fluor_suite2P’- looking at the all cells fluorescence and stained cells fluorescence after visually verify that cells are positively stained. For behavior.
‘run_ Zstack_get_fluor_suite2P’- looking at Z-stack fluorescence data, and compare non-stained to stained. 
'find_stained_cells' used to identified the stained cells

*** More detailed steps:
1)	Run ‘get_data_suite2P’. make sure to define mice and file locations. 
Set all initial save options to 1. 
Obeserve the contour with 647 label figure ‘F647-1’ with imageJ. It might be that manual correction to suite2p positive cells is needed, to look at cells that are positive to staining (‘file->load processed data., find the relevant ‘suite2p’ folder and open the ‘stat.npy’ file) 
If so- edit with suite2p. make sure to save to *.mat (file->save to mat file)
Re-run ‘get_dat_suite2P’ with save_contour=1 new_contour=1 
After this is done, add the identify cell to the array: ’stained_cell_ind’
To compare few sessions- few iterations are needed to find the same cells in all sessions. 
I recommend save the combined tif from imageJ for future needs



Figure S1: 

FP with reward conditioning

Plot_RewardConditioning_processed


