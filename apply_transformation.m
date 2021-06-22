function apply_transformation
% load the tform and apply it 
clear FIT
loadfit=0;
% load fixed figure
%fixedimage=imread(['C:\Users\anatk\Documents\Light_sectioning\WT316RR_LGS\885um\885um_registered_invivo2GRIN.tif']);
%fixedimage=imread(['C:\Users\anatk\Documents\Light_sectioning\WT_242_LGS\740um\740um_registered_invivo2GRIN.tif']);
%fixedimage=imread('C:\Users\anatk\Documents\Light_sectioning\WT36R_LGS\485um\485um_registered_invivo2GRIN.tif');

fixedimage=imread('C:\Users\anatk\Documents\Light_sectioning\WT35L_LGS\790um\790um_registered_invivo2GRIN_NR.tif');

%fixedimage=imread('C:\Users\anatk\Documents\Light_sectioning\WT58N_LGS\655um\655um_registered_invivo2GRIN.tif');
%fixedimage=imread('C:\Users\anatk\Documents\Light_sectioning\WT2170\1290um\section_12_Invivo2GRIN940nm.tif');

% load moving figure
%movingdimage=imread(['C:\Users\anatk\Documents\Light_sectioning\WT316RR_LGS\316RR_invivo\MAX\MAX_VIPWT_316RR_94mW_500_256_13_Zstack_invivo-885um.tif']);
%movingimage=imread(['C:\Users\anatk\Documents\Light_sectioning\WT_242_LGS\242_invivo\VIP242_Max\MAX_VIP242_87uW_740um_081618.tif']);
%movingimage=imread('C:\Users\anatk\Documents\Light_sectioning\WT36R_LGS\WT36R_invivo\WT36R_invivo5\WT36R_MC_MAX5\section_1_-485nm_v2_nrMC_max.tif');
movingimage=imread('C:\Users\anatk\Documents\Light_sectioning\WT35L_LGS\WT35L_invivo1\WT35L_MC_MAX\section_17_-790nm_v2_nrMC_max.tif');

%movingimage=imread('C:\Users\anatk\Documents\Light_sectioning\WT58N_LGS\58N_invivo\WT58N_MC_MAX\s-655nm_v2_nrMC_max.tif');
%movingimage=imread('C:\Users\anatk\Documents\Light_sectioning\WT2170\WT2170_invivo\WT2170_invivo4\WT2170_maxMC4\section_11_-1275nm_v2_nrMC_max.tif');


% load transform 
switch loadfit
    case 1
        %FIT_path='C:\Users\anatk\Documents\Light_sectioning\WT316RR_LGS\fit_files\316RR_invivo1p1_885_FIT_NRS.mat';
        FIT_path='C:\Users\anatk\Documents\Light_sectioning\WT_242_LGS\fit_files\242_invivo1p1_740_FIT_NRS.mat';
        
        load(FIT_path);
        Fix_tform=FIT.mytform;
        % apply transformation
        transformed_image=uint16(imwarp(movingimage,Fix_tform,'OutputView',imref2d(size(fixedimage))));
        
    case 0
        %R=flip(imrotate(fliplr(movingimage),90));
      %  R=fliplr(50*movingimage);
      R=3*movingimage;
        [movingPoints,fixedPoints]= cpselect(R, fixedimage,'Wait',true);
        mytform = fitgeotrans(movingPoints, fixedPoints, 'NonreflectiveSimilarity');
         transformed_image=uint16(imwarp(R,mytform,'OutputView',imref2d(size(fixedimage))));
         [MOVINGREG.DisplacementField,transformed_image2] = imregdemons(transformed_image,fixedimage,100,'AccumulatedFieldSmoothing',2,'PyramidLevels',1);
   end

figure
imshowpair(transformed_image,fixedimage,'montage')

%% save
%cd ('C:\Users\anatk\Documents\Light_sectioning\WT316RR_LGS\885um')
%cd ('C:\Users\anatk\Documents\Light_sectioning\WT_242_LGS\740um')
%cd ('C:\Users\anatk\Documents\Light_sectioning\WT36R_LGS\485um')
cd ('C:\Users\anatk\Documents\Light_sectioning\WT35L_LGS\790um')
%cd ('C:\Users\anatk\Documents\Light_sectioning\WT58N_LGS\655um')
%cd('C:\Users\anatk\Documents\Light_sectioning\WT2170\1290um')

imwrite(transformed_image,['790um_MAX_invivo2GRIN.tif']);
imwrite(transformed_image2,['790um_MAX_invivo2GRIN_NR.tif']);
