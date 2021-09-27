%% after running suite2P, a Fall.mat file is generated 
% Anat Kahan Cell Reports 2021
% initialize variables 
stained_cell_ind=[];
clear M cell_ind mytform Yregistered Y
read_YMat=0;% set to 1 at the first time running the script- save the in vivo data. then it can go back to zero
new_contours=0;
new_647_contours=0;
new_555_contours=1;
save_contours=1;
new_647=0;
new_555=1;
STAINED=3; % 1 is for 647. 2 is for non-stained. 3 is for 555
% set path
path='C:\Users\anatk\Documents\Light_sectioning\';
 A=[];    
     
%mouse='WT36R_LGS'; exp='Zstack'
mouse='Str39_LGS' ; exp='behavior_com'

switch mouse
    case 'WT36R_LGS'
        M_num='36R';
        % sess=755; fixed_647_section=26; fixed_GCaMP_section=87;stained_cell_ind=16;
         sess=485; fixed_647_section=40; fixed_GCaMP_section=49;stained_cell_ind=[3 4 21];
       % sess=395; fixed_647_section=61; fixed_GCaMP_section=24;stained_cell_ind=[18 1 7 50];%[4,14,19,20]- sparse mode;
        
        fixed_path=['WT36R_LGS\' num2str(sess) 'um\'];
        exp='Zstack'; Z=sess;
    case 'Str39_LGS'
        M_num='Str39';
        
        sess=1; fixed_647_section=28; fixed_GCaMP_section=48;%stained_cell_ind=[];%[4,14,19,20]- sparse mode;
        
        fixed_path=['Str39_LGS\' num2str(sess) 'um\'];
         Z=870;  
     

         stained_cell_ind=[3 9 36 39 29 82 87]; 
         stained_cell_ind_555=[82 87]; 

         non_stained_cell_ind=[38 68 12 69 17 80 18 ];
         A=intersect(stained_cell_ind,non_stained_cell_ind);
         if ~isempty(A)
             disp('non stained and stained has a common cell- CHECK!!!!!')
         end
         
         switch STAINED
            case 1;  stained_cell_ind=stained_cell_ind(~isnan(stained_cell_ind));
            case 2;  stained_cell_ind=non_stained_cell_ind(~isnan(non_stained_cell_ind));
            case 3;  stained_cell_ind=stained_cell_ind_555(~isnan(stained_cell_ind_555));
          
        end

         
end

switch exp
    case 'behavior'
        this_dir='all_tiff_to_align';
    case 'behavior_com'
        this_dir='all_tiff_to_align_combined';
        
    case 'Zstack'
        this_dir='tiff_Zstack_to_align';
end

cd ([path mouse '\' this_dir '\sess' num2str(sess) '\suite2p\plane0\'])


switch mouse
    case 'WT36R_LGS'
        fixedimage=imread([path mouse '\' num2str(Z) 'um\' num2str(Z) 'um_registered_invivo2GRIN.tif']);
    case 'Str39_LGS'
        fixedimage=imread([path mouse '\' num2str(Z) 'um\' num2str(Z) 'um_registered_invivo2GRIN.tif']);

end

M=load('Fall.mat');
cell_ind=find(M.iscell(:,1)>0);


%% plot cells pixels
IM=10*uint16(M.ops.meanImg);
figure
%im2=imshow(IM);
subplot(1,2,1);imshow(IM);title(['behavior sess ' num2str(sess)]);subplot(1,2,2);imshow(fixedimage);title('Z-stack')

% get info from Fall Matrix - the matlab version of the data created by
% suite2P, get the identified cells pixels (xpix, ypix) and creates contour map
nim=creates_contours_figure(cell_ind, M,mouse,sess);

% creates a figure of the stained positive cells, if defined 
if ~isempty(stained_cell_ind)
    Stained_nim=creates_contours_figure(stained_cell_ind', M,mouse,sess);
end

% now fit the image of the behavior to the Z-stack (fixedimage)
fullpath=[path mouse '\' this_dir '\sess' num2str(sess) '\suite2p\plane0\FIT.mat' ];
clear mytform movingPoints fixedPoints
if ~exist(fullpath)
    % define similarity points    
    [MovingPoints,FixedPoints] = cpselect(IM, fixedimage,'Wait',true);
    
    %% transformation. choose the type of transformation carefully
    % movingPoints=evalin('base','movingPoints');
    % fixedPoints=evalin('base','fixedPoints');
    mytform = fitgeotrans(MovingPoints, FixedPoints, 'NonreflectiveSimilarity');
    
    save(fullpath,'mytform')
else
    load(fullpath) % loads the 'mytform'
end

%% apply the tranformation to the moving figure
Jregistered = imwarp(IM,mytform,'OutputView',imref2d(size(fixedimage)));
figure; imshowpair(Jregistered,fixedimage); title([mouse ' sess '  num2str(sess)])
%  figure; subplot(1,2,1);imshow(Jregistered);subplot(1,2,2);imshow(fixedimage);
%% apply transformation to the contour image
Cregistered = imwarp(nim,mytform,'OutputView',imref2d(size(fixedimage)));
if ~isempty(stained_cell_ind)
   SCregistered = imwarp(Stained_nim,mytform,'OutputView',imref2d(size(fixedimage)));   
end


cd('../..')

% apply the tranformation to the contour figure 
imwrite(uint16(Cregistered),'contours.tif','tiff');
imwrite(uint16(Jregistered),['in_vivo_MIPs_sess' num2str(sess) '.tif'],'tiff');
figure; imshowpair(Cregistered,fixedimage);title([mouse ' sess '  num2str(sess)])
if ~isempty(stained_cell_ind)
    if save_contours
        switch STAINED
            case 1; imwrite(uint16(SCregistered),['Stained_contours_sess' num2str(sess) '.tif'],'tiff');
            case 2; imwrite(uint16(SCregistered),['non_Stained_contours_sess' num2str(sess) '.tif'],'tiff');
            case 3; imwrite(uint16(SCregistered),['Stained_555_contours_sess' num2str(sess) '.tif'],'tiff');
  
        end
    end
  figure; imshowpair(SCregistered,fixedimage)  
   imwrite(fixedimage,'fixed_image_mean.tif','tiff');
  title([mouse ' sess '  num2str(sess)])
end
% save the contour X times of frames
if new_contours
    mulitple_save_image_contours_to_behavior(uint16(Cregistered),mouse,sess,size(M.F,2),exp,'contours')
end
if new_647_contours || new_555_contours
    if ~isempty(stained_cell_ind)
        switch STAINED
            case 1 ; mulitple_save_image_contours_to_behavior(uint16(SCregistered),mouse,sess,size(M.F,2),exp,'stained_contours')
            case 2 ; mulitple_save_image_contours_to_behavior(uint16(SCregistered),mouse,sess,size(M.F,2),exp,'non_stained_contours') 
             case 3 ; mulitple_save_image_contours_to_behavior(uint16(SCregistered),mouse,sess,size(M.F,2),exp,'stained_555_contours')
 
        end
    end
end


%%apply the transformation on sess invivo data
switch exp
    case 'behavior'
        Y=read_file(['uint16YMat_' num2str(sess-1) '.tif']);
    case 'behavior_com'
        %Y=read_file(['uint16YMat_' num2str(sess-1) '.tif']);
        Y=read_file(['s-' num2str(sess) '.tif']);
        Y2=read_file(['s-' num2str(sess+1) '_transformed.tif']);
    case 'Zstack'
        Y=read_file(['uint16YMat.tiff']);
end

% first time running the script- save the in vivo data 
if read_YMat
    if ~exist('YMat_registered','dir'); mkdir ('YMat_registered');end
    cd('YMat_registered')
    for i=1:size(Y,3)
        k=i;
        switch exp
            case 'behavior_com'
              %  Y=Y1; 
               Y=Y2; k=2391+i;
        end
        YR=imwarp(Y(:,:,i),mytform,'OutputView',imref2d(size(fixedimage)));
        %    otherwise
        %        YR=Y(:,:,i);
        %end
        Yregistered(:,:,i) =YR;
       
        imwrite(YR,['Y' num2str(k) '.tif']);
    end
end
%cd('..')
% now contour and TMat_register folder in be improrted as 'squence' image in imageJ, and merged using 'color - merge'

 %% read the fixed tissue file. for 2170 the best imaging was from confocal 
 
 switch mouse
      case {'WT36R_LGS'}
         fix_GCaMP=imread([path fixed_path 'section_' num2str(fixed_GCaMP_section) '.tif']);
        % fix_647=imread([path fixed_path 'section_' num2str(fixed_647_section) '.tif']);    
            fix_647=imread([path mouse  '\WT36R_fixed_2P\1100_15mW_adjusted_042821\1100_15mW_adjusted_042821B46.tif']);  

         fullpath=[path mouse '\fit_files\' M_num '_processed_940nm_' num2str(Z) '_FIT_NRS.mat' ];

       case {'Str39_LGS'}
         fix_GCaMP=imread([path mouse '\' num2str(Z) 'um\section_' num2str(fixed_GCaMP_section) '_GRIN940nm.tif']);
         fix_647=imread([path mouse '\confocal_processed\AVG_C3.tif']); % 
         fix_555=imread([path mouse '\confocal_processed\AVG_C2.tif']); % 
         fullpath=[path mouse '\fit_files\' M_num(end-1:end) '_Confocal2twoP_' num2str(Z) '_FIT_NRS.mat' ];

 end

if size(fix_GCaMP,3)>1; fix_GCaMP=fix_GCaMP(:,:,2); end
if size(fix_647,3)>1; fix_647=fix_647(:,:,2); end

% load the tform for the fixed tissue and staining
load(fullpath);
Fix_tform=FIT.mytform;
F647R=uint16(imwarp(fix_647,Fix_tform,'OutputView',imref2d(size(fixedimage))));
switch mouse
    case 'Str39_LGS'
        F555R=uint16(imwarp(fix_555,Fix_tform,'OutputView',imref2d(size(fixedimage))));
end
FGcaMPR=uint16(imwarp(fix_GCaMP,Fix_tform,'OutputView',imref2d(size(fixedimage))));

%figure; imshow(fix_647)

figure
subplot(1,2,1);imshowpair(Cregistered,F647R);
%% figure to show cell contour with invivo 
% get contours using nim processing
subplot(1,2,2); imshowpair(Cregistered,Jregistered)
title (['sess ' num2str(sess)])
% save 647 to folder X number of frames

cd ([path mouse '\' this_dir '\sess' num2str(sess) '\'])
imwrite(FGcaMPR,'FGCaMP_1.tif');

if new_647
    if ~exist('647_registered','dir'); mkdir ('647_registered');end
    imwrite(F647R,'F647_1.tif');
    
    cd('647_registered')
    for i=1:size(Y,3)*2
        imwrite(F647R,['F647_' num2str(i) '.tif']);
    end
    cd('..')
end

if new_555
    if ~exist('555_registered','dir'); mkdir ('555_registered');end
    imwrite(F555R,'F555_1.tif');
    
    cd('555_registered')
    for i=1:size(Y,3)*2
        imwrite(F555R,['F555_' num2str(i) '.tif']);
    end
    cd('..')
end

if ~isempty(A)
    disp('non stained and stained has a common cell- CHECK!!!!!')
end
1