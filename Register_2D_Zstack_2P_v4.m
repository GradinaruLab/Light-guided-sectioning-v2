function Register_2D_Zstack_2P_v4 
% load 1 plane data, and fit by manually finding the same neurons
% v4 includes a non rigid correction to the affine transformation 

close all
clear AfixedImage registered Jregistered fixedImage movingImage AmovingImage

%1 is in vivo vs processed GRIN; GRIN processed is fixed;
%2 is procesessd, GRIN vs tissue side; GRIN processed is fixed; 
%3 is confocal to 2P, tissue side,, after  tissue was registered 
comp=2
i=5; % index of plane 

file_name='LGS registration table home'; % 3 is with estrus state refered to proestrus
path='D:\Data_Glab\Light_sectioning\';
path='C:\Users\anatk\Documents\Light_sectioning\';
full_path=[path file_name '.xlsx'];
%[NUM,TXT,RAW]=xlsread(full_path,'VS AK');
[NUM,TXT,RAW]=xlsread(full_path);

%mouse='WT36R_LGS';num_reg_points=5;invivo='WT36invivo1';this_mouse_ind=find(contains(TXT(:,1),invivo));
mouse='WT36R_LGS';num_reg_points=7;
%mouse='Str39_LGS';num_reg_points=4;

this_mouse_ind=find(contains(TXT(:,1),mouse(1:end-4)));

% extract parameters from file
all_Z=-1*NUM(this_mouse_ind+4:this_mouse_ind+4+num_reg_points-1,2);
all_s2=NUM(this_mouse_ind+4:this_mouse_ind+4+num_reg_points-1,3);
all_s1=NUM(this_mouse_ind+4:this_mouse_ind+4+num_reg_points-1,5);
all_s1100=NUM(this_mouse_ind+4:this_mouse_ind+4+num_reg_points-1,7);


inds= strfind(mouse,'_');
s_ind=strfind(mouse,'Str');
if length(inds)>1 
    M_num=mouse(inds(1)+1:inds(2)-1);
    mouse_title=[mouse(1:inds(1)-1) ' ' M_num ' ' mouse(inds(2)+1:end)];
elseif ~isempty(s_ind)
    M_num=mouse(4:inds(1)-1);
    mouse_title=[mouse(1:3) ' ' M_num ];
else
    M_num=mouse(3:inds(1)-1);
    mouse_title=[mouse(1:2) ' ' M_num ' ' mouse(inds(1)+1:end)];
end

%% those parameters are taken from 'LGS registration table'

Z=all_Z(i); % IN-VIVO
choose_section2=all_s2(i); % PROCESSED GRIN side
choose_section3=all_s1(i); % PROCESSED Tissue side
choose_section_T1100=all_s1100(i);
% registration method:
method_fixed='NRS';

if ~isempty(strmatch(mouse,'Str39_LGS')) && Z==870 ; Z2=600; else Z2=Z; end
% register confocal directly to invivo 
%if ~isempty(strmatch(mouse,'Str39_LGS')) && Z==870 ; choose_section2=choose_section_T1100; end

switch comp
    case 1
        fullpath=[path mouse '\fit_files\' M_num '_invivo1_' num2str(Z) '_FIT_NRS.mat'];
    case 2
        fullpath=[path mouse '\fit_files\' M_num '_processed_940nm_' num2str(Z) '_FIT_' method_fixed '.mat' ];
    case 3
        fullpath=[path mouse '\fit_files\' M_num '_Confocal2twoP_' num2str(Z) '_FIT_' method_fixed '.mat' ];
end

if exist(fullpath)
   load(fullpath); 
end

%     

clear atpath
%% fixed, through GRIN lens path, 940nm:
atpath{1}=[TXT{this_mouse_ind+2,4} '\'];
%% in vivo path: pay attaention here that 36R has 'combined' and 'combined not adjusted' folders
if strfind(TXT{this_mouse_ind+2,2},'MAX')
    atpath{2}=[TXT{this_mouse_ind+2,2} '\s-' num2str(Z2) 'um_v2_nrMC_max.tif'];
else
    atpath{2}=[TXT{this_mouse_ind+2,2} '\s-' num2str(Z2) 'um_v2_nrMC_mean.tif'];
end
%% fixed, through Tissue path, 940nm:
atpath{3}=[TXT{this_mouse_ind+2,6} '\'];
%% fixed, through Tissue path, 1100nm:
atpath{4}=[TXT{this_mouse_ind+2,8} '\'];

atpath{6}=atpath{1};%% might need adjusted for 36R: %atpath{6}=['D:\DATA_Glab\Light_sectioning\WT36R_LGS\WT36R_fixed_2P\WT36R_GRIN_1100nm_350_512_125_27mW\adjusted\'];

switch comp
    case 1 % in vivo vs processed. In vivo is the fixed
        % load processed/Fixed GRIN imaging data
        disp('load ''processed GRIN lens'' data')
        % 02/21/2019 I changed moving/fixed image
          % 01/31/21- I changed back moving/fixed- processed GRIN should be
          % the reference (when available)
        fixedImage=imread([atpath{1} 'section_' num2str(choose_section2) '.tif']);
        switch mouse
            case 'Str39_LGS'
                fixedImage=imread([atpath{4} 'AVG_C1.tif']);
        end
                 
        % load in vivo data
        disp('load ''in vivo'' data')
        tmpImage=imread(atpath{2});
        movingImage=imrotate(tmpImage,0);

    case 2 % procesessd, GRIN vs tissue side,940nm   

        % load processed GRIN imaging data
        disp('load ''processed GRIN lens'' data')
        fixedImage=imread([atpath{1} 'section_' num2str(choose_section2) '.tif']);
        
        if size(fixedImage,3)>1; fixedImage=fixedImage(:,:,2);end
  
        % load processed tissue imaging data
        disp('load ''processed Tissue lens'' data')
        
        if strcmp(mouse,'Str39_LGS')
                %movingImage=imread([atpath{4} 'AVG_C1.tif']);
 
         movingImagetmp=imread([atpath{3} 'section_' num2str(choose_section3) '.tif']);
           movingImage=flip(movingImagetmp);
        else
            movingImagetmp=imread([atpath{3} 'section_' num2str(choose_section3) '.tif']);
            movingImage=fliplr(movingImagetmp);
        end
        
        if size(movingImage,3)>1; movingImage=movingImage(:,:,2);end
        
      case 3 % procesessd, tissue side, 2P vs confocal

        % load processed GRIN imaging data
        disp('load ''processed Tissue, 2 photon'' data')
        switch mouse
            case 'Str39_LGS'
                fixedImage=imread([path mouse '\' num2str(Z) 'um\' num2str(Z) 'um_registered_Tissue2GRIN.tif']);
            otherwise
                fixedImage=imread([atpath{3} 'section_adjusted_' num2str(choose_section3) '.tif']);
        end
        
        % load processed tissue imaging data
        disp('load ''processed Tissue, confocal'' data')
        switch mouse
            case 'WT36R_LGS'
                movingImagetmp=imread([atpath{4} 'WT36R_Tissue_0to320_FOXP2_647_GFP_488_Zstack_10X_LGS.tiff_files\488\WT36R_Tissue_0to320_FOXP2_647_GFP_488_Zstack_10X_LGS_b0v0t0z' num2str(choose_section4-1) 'c0x0-2048y0-2048.tiff']);
                movingImagetmp2=fliplr(movingImagetmp);
                movingImage=imrotate(movingImagetmp2,90);
                movingImage=im2uint16(movingImage);% increases image quality
                movingImage=movingImage(:,:,2);
            case 'Str39_LGS'
            %    movingImage=imread([path mouse '\' num2str(Z) 'um\C1-AVG_Str39_Airyscan_488_cfos555_ARC647_25X_dip_SeeDB_Zstack_B_Airyscan Processing.tif']);
              movingImage=imread([path mouse '\confocal_processed\AVG_C1.tif']);
                
        end
end

%% brightness adjustment 
BrF(1)=20;
BrF(2)=8;
[AmovingImage, AfixedImage]=adjust_brightness(uint16(movingImage),fixedImage,BrF(1),BrF(2),1);

if ~exist(fullpath)
    %% define similarity points
   [movingPoints,fixedPoints]= cpselect(AmovingImage, AfixedImage,'Wait',true);

    %% transformation. choose the type of transformation carefully
%     movingPoints=evalin('base','movingPoints');
%     fixedPoints=evalin('base','fixedPoints');
   mytform = fitgeotrans(movingPoints, fixedPoints, 'NonreflectiveSimilarity');
    % mytform = fitgeotrans(movingPoints, fixedPoints, 'lwm');
else
    mytform=FIT.mytform;
end

%% apply the tranformation to the moving figure 
Jregistered = imwarp(AmovingImage,mytform,'OutputView',imref2d(size(AfixedImage))); 
Jregistered=uint16(Jregistered);

[MOVINGREG.DisplacementField,Jregistered2] = imregdemons(Jregistered,AfixedImage,100,'AccumulatedFieldSmoothing',2,'PyramidLevels',1);

% looking for the new tform
[optimizer, metric] = imregconfig('multimodal');

%Affine transformation consisting of translation, rotation, scale, and shear
NR_tform = imregtform(Jregistered2,Jregistered,'affine',optimizer,metric);
Jregistered3 = imwarp(Jregistered,NR_tform,'OutputView',imref2d(size(Jregistered))); 

ssim(AfixedImage,Jregistered)
ssim(AfixedImage,uint16(Jregistered2))
ssim(Jregistered3,Jregistered2)

figure
subplot(1,2,1)
imshowpair(Jregistered3,Jregistered2,'montage')
xlabel('manual')
subplot(1,2,2)
imshowpair(AfixedImage,Jregistered2,'montage')
xlabel('manual with NR')

fake1=zeros(size(AfixedImage,1));

figure
subplot (1,4,3), imshowpair(AfixedImage,fake1),title('adj.B fixed Image ')
subplot (1,4,1), imshow(AmovingImage), title('adj.B moving Image ');set(gca,'xtick',[0 175 350]);set(gca,'xticklabel',{'0','175','350'});
subplot (1,4,4), imshowpair(AfixedImage,Jregistered2), title('fixed and registered Image')
subplot (1,4,2), imshowpair(fake1,Jregistered2), title('registered moving Image')

[FIT.this_angle, FIT.this_scale]= LGS_get_angle_and_scale(mytform);
FIT.mytform=mytform;
FIT.BrF=BrF;
FIT.NR_tform=NR_tform;

switch comp
    case {1, 1.1}
        text(-250,600,[  'section ' num2str(choose_section2) ', ' num2str(Z) 'um ' mouse_title ', NRS'])
    case {2,2.1}
        text(-350,600,[ 'section ' num2str(choose_section2) ',to section' num2str(choose_section3) ', ' mouse_title, ', ' method_fixed])
    case 3
        text(-350,600,[ num2str(Z) 'um ' mouse_title, ', ' method_fixed])
end

figure, imshowpair(AfixedImage,Jregistered2), title('fixed image and registered Image')

switch comp
    case 1
        save([path mouse '\fit_files\' M_num '_invivo1_' num2str(Z) '_FIT_NRS'],'FIT')
        sind=strfind(atpath{2},'\');
        imwrite(Jregistered,[atpath{2}(1:sind(end))   num2str(Z)  'um_registered_invivo2GRIN.tif'])
        imwrite(Jregistered2,[atpath{2}(1:sind(end))   num2str(Z)  'um_registered_invivo2GRIN_NR.tif'])
  
          if ~exist([path mouse '\'   num2str(Z) 'um\']); mkdir([path mouse '\'   num2str(Z) 'um\']); end
          cd ([path mouse '\'   num2str(Z) 'um\']); 
          imwrite(Jregistered,[num2str(Z)  'um_registered_GRIN2invivo.tif'])
          imwrite(Jregistered2,[num2str(Z)  'um_registered_GRIN2invivo_NR.tif'])
    case 2
        save([path mouse '\fit_files\' M_num '_processed_940nm_' num2str(Z) '_FIT_' method_fixed ],'FIT')
        imwrite(Jregistered,[atpath{3} 'section_registered_' num2str(choose_section3) '.tif'])
          imwrite(Jregistered2,[atpath{3} 'section_registered_' num2str(choose_section3) '_NR.tif'])
     
        if ~exist([path mouse '\'   num2str(Z) 'um\']); mkdir([path mouse '\'   num2str(Z) 'um\']); end
        imwrite(Jregistered,[path mouse '\'   num2str(Z) 'um\' num2str(Z) 'um_registered_Tissue2GRIN.tif'])
        imwrite(Jregistered2,[path mouse '\'   num2str(Z) 'um\' num2str(Z) 'um_registered_Tissue2GRIN_NR.tif'])
      
        imwrite(AfixedImage,[path mouse '\'   num2str(Z) 'um\' 'section_' num2str(choose_section2) '_GRIN940nm.tif'])
    case 3 % the matrix saved is in relate to **adjusted** tissue, 2P!!!!  
        save([path mouse '\fit_files\' M_num '_Confocal2twoP_' num2str(Z) '_FIT_' method_fixed ],'FIT')
        imwrite(Jregistered,[path mouse '\'   num2str(Z) 'um\' num2str(Z) 'um_registered_Confocal2twoP.tif'])
        imwrite(Jregistered2,[path mouse '\'    num2str(Z) 'um\' num2str(Z) 'um_registered_Confocal2twoP_NR.tif'])

end

