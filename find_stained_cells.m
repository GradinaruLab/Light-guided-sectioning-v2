function find_stained_cells
% define the positive stained cells for LiGS


% read the contours, after running 'get_data_suite2P'
%cd ('C:\Users\anatk\Documents\Light_sectioning\Str39_LGS\all_tiff_to_align_combined\sess1')
cd ('C:\Users\anatk\Documents\Light_sectioning\WT36R_LGS\tiff_Zstack_to_align\sess485')

C=imread('contours.tif');

% % 
% M=imread('870um_ARC_crop.tif');
% C=imread('870um_contours_crop.tif');

% read the stained image, after using
% 'apply_transformation_on_stained_image'
%cd ('C:\Users\anatk\Documents\Light_sectioning\Str39_LGS\870um')
cd ('C:\Users\anatk\Documents\Light_sectioning\WT36R_LGS\755um')
M=imread('F647_755um_C.tif');
F=12;
%M=imread('F555_870um_C.tif');

M=10*M;
M2=M>F*std2(M);
figure
imshow(M2)

figure
imshowpair(M2,C)



