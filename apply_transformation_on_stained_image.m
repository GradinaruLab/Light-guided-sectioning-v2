function apply_transformation_on_stained_image
% load the tform for the fixed tissue and staining
file_name='LGS registration table home'; % 3 is with estrus state refered to proestrus
%path='D:\Data_Glab\Light_sectioning\';
path='C:\Users\anatk\Documents\Light_sectioning\';
clear FIT

full_path=[path file_name '.xlsx'];
[NUM,TXT,RAW]=xlsread(full_path);

comp=2 %1 is in vivo vs processed; 1.1 is invivo vs processed, processed is fixed; 2 is procesessd, GRIN vs tissue side; ; 2.1 is procesessd, GRIN vs tissue side, 1100nm;3 is confocal to 2P, tissue side; 4 is invivo to 2P
i=3;

% brains from 05-06/2019
%mouse='WT58N_LGS'; num_reg_points=4;this_mouse_ind=find(contains(TXT(:,1),mouse(1:end-4)));
%mouse='Drd1_1N';num_reg_points=4; this_mouse_ind=find(contains(TXT(:,1),mouse));
%mouse='WT35L_LGS'; num_reg_points=4;this_mouse_ind=find(contains(TXT(:,1),mouse(1:end-4)));
%mouse='WT36R_LGS';num_reg_points=5;invivo='WT36invivo1';this_mouse_ind=find(contains(TXT(:,1),invivo));
mouse='WT36R_LGS';M_num='36R' ; num_reg_points=6;this_mouse_ind=find(contains(TXT(:,1),mouse(1:end-4)));
%mouse='WT316RR_LGS';num_reg_points=6;this_mouse=strfind(TXT(:,1),mouse(1:end-4));
%mouse='Drd1a_1816L_LGS';num_reg_points=4;this_mouse_ind=find(contains(TXT(:,1),mouse(1:end-4)));
%mouse='WT_242_LGS';num_reg_points=6;this_mouse_ind=find(contains(TXT(:,1),mouse(1:end-4)));
%mouse='WT2170_LGS';M_num=2170;num_reg_points=4;
%mouse='WT2174_LGS';M_num=2174;num_reg_points=4; 
%mouse='Str39_LGS';M_num=39;num_reg_points=4; ch=[1 2 3];
ch=[3]; 
%fixed_path='WT2170_processed_confocal\WT2170_GC488_ARC_647_10X_zoom2p5_3um_step_10X.tiff_files\';

this_mouse_ind=find(contains(TXT(:,1),mouse(1:end-4)));
fixed_path=TXT{this_mouse_ind+2,8};
all_Z=-1*NUM(this_mouse_ind+4:this_mouse_ind+4+num_reg_points-1,2);
all_647_section=NUM(this_mouse_ind+4:this_mouse_ind+4+num_reg_points-1,7);
Z=all_Z(i); % IN-VIVO
fixed_section=all_647_section(i);
%
switch comp
    case 2
          fixedimage=imread([path mouse '\' num2str(Z) 'um\'  num2str(Z) 'um_registered_Tissue2GRIN.tif']);
        %fixedimage=imread([path mouse(1:end-4) '\' num2str(Z) 'um\section_' num2str(fixed_section) '_Tissue2GRIN940nm.tif']);

        fullpath=[path mouse '\fit_files\' M_num '_processed_940nm_' num2str(Z) '_FIT_NRS.mat' ];
    case 3
        fixedimage=imread([path mouse '\' num2str(Z) 'um\section_' num2str(fixed_section) '_Confocal2TwoP.tif']);

        fullpath=[path mouse '\fit_files\' num2str(M_num) '_Confocal2twoP_' num2str(Z) '_FIT_NRS.mat' ];        
        %       fullpath=[path mouse '\fit_files\' M_num '_Confocal2twoP_' num2str(Z) '_FIT_NRS.mat' ];

end

% load 647 figure
%fix_647=imread([fixed_path '\section_' num2str(fixed_section) '.tiff']);
switch mouse
    case 'Str39_LGS'
        fix_647=imread([path mouse  '\confocal_processed\AVG_C3.tif']);
        if find(ch==2)
            fix_555=imread([path mouse  '\confocal_processed\AVG_C2.tif']);
        end
        if find(ch==1)
            fix_488=imread([path mouse  '\confocal_processed\AVG_C1.tif']);
        end
        if find(ch==2)
            if size(fix_555,3)>1; fix_555=fix_555(:,:,2); end % needed if images were taken with Zeiss confocal
        end
    case 'WT36R_LGS'
      if i==5 
       fix_647=imread([path mouse  '\WT36R_fixed_2P\1100_15mW_adjusted_042821\1100_15mW_adjusted_042821B46.tif']);  
      end
      if i==3
        fix_647=imread([path mouse  '\WT36R_fixed_2P\1100_15mW_adjusted_042821\1100_15mW_adjusted_042821B30.tif']);  
      end
end

if size(fix_647,3)>1; fix_647=fix_647(:,:,2); end % needed if images were taken with Zeiss confocal

figure; imshow(10*uint16(fix_647))

% load transform 
load(fullpath); 
Fix_tform=FIT.mytform;

% apply transformation
F647R=uint16(imwarp(fix_647,Fix_tform,'OutputView',imref2d(size(fixedimage))));
if find(ch==2)
    F555R=uint16(imwarp(fix_555,Fix_tform,'OutputView',imref2d(size(fixedimage))));
end
if find(ch==1)
    F488R=uint16(imwarp(fix_488,Fix_tform,'OutputView',imref2d(size(fixedimage))));
end
%% save
imwrite(F647R,[path mouse '\' num2str(Z) 'um\F647_' num2str(Z) 'um_C.tif']);
if find(ch==2)
    imwrite(F555R,[path mouse '\' num2str(Z) 'um\F555_' num2str(Z) 'um_C.tif']);
end
if find(ch==1)
    imwrite(F488R,[path mouse '\' num2str(Z) 'um\F488_' num2str(Z) 'um_C.tif']);
end