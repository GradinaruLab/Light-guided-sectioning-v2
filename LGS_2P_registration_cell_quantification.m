function [image_score,cell_score,is_cell,shuffeled_score]=LGS_2P_registration_cell_quantification(mouse,i,thresh)
% following Register_2D_Zstack_2P_v3. check how many cells can be found 

%close all 
clear AfixedImage registered Jregistered fixedImage movingImage AmovingImage atpath 
clear image_score

comp=1.1 ;%1.1 is invivo vs processed fixed; 2.1 is procesessd, GRIN vs tissue side, 1100nm;




% brains from 05-06/2019
%%%%%%%% read Z and indexes from file
file_name='LGS registration table home'; % 3 is with estrus state refered to proestrus

path='C:\Users\anatk\Documents\Light_sectioning\';
full_path=[path file_name '.xlsx'];
%[NUM,TXT,RAW]=xlsread(full_path,'VS AK');
[NUM,TXT,RAW]=xlsread(full_path);

% brains from 05-06/2019
invivo=[];
if nargin == 0
  % mouse='WT58N_LGS'; 
    %mouse='Drd1_1N'; num_reg_points=4;this_mouse_ind=find(contains(TXT(:,1),mouse));
    %mouse='WT35L_LGS'; num_reg_points=5;this_mouse_ind=find(contains(TXT(:,1),mouse(1:end-4)));
   %mouse='WT36R_LGS';
    
   % mouse='WT36R_LGS';
    mouse='WT316RR_LGS';
    %mouse='Drd1a_1816L_LGS';num_reg_points=4;this_mouse_ind=find(contains(TXT(:,1),mouse(1:end-4)));
    %mouse='WT_242_LGS';num_reg_points=6;this_mouse_ind=find(contains(TXT(:,1),mouse(1:end-4)));  
    %mouse='WT2170_LGS';num_reg_points=3;this_mouse_ind=find(contains(TXT(:,1),mouse(1:end-4)));
    thresh=0.5;
    i=5; % index of plane
end
this_mouse_ind=find(contains(TXT(:,1),mouse(1:end-4)));
if i<3
    radius=45;
else
    radius=30;
end
switch mouse
    case 'WT58N_LGS'; num_reg_points=5; radius=30;
    case 'Drd1_1N'; num_reg_points=4;
    case 'WT35L_LGS'; num_reg_points=5;
    case 'WT36R_LGS';num_reg_points=5;
    case 'WT316RR_LGS';num_reg_points=6;if i<3 radius=55; end; if i==3 || i==4 radius=45; end
    case 'Drd1a_1816L_LGS';num_reg_points=4;
    case 'WT_242_LGS';num_reg_points=6;
    case 'WT2170_LGS';num_reg_points=3;
end



all_Z=-1*NUM(this_mouse_ind+4:this_mouse_ind+4+num_reg_points-1,2);
all_s2=NUM(this_mouse_ind+4:this_mouse_ind+4+num_reg_points-1,3);   
Z=all_Z(i); % IN-VIVO
choose_section2=all_s2(i); % PROCESSED GRIN side
tmp_invivo_registered=imread([path mouse '\'   num2str(Z) 'um\' num2str(Z) 'um_registered_invivo2GRIN.tif']);
fixed_GRIN=imread([path mouse '\'   num2str(Z) 'um\' 'section_' num2str(choose_section2) '_GRIN940nm.tif']);
[MOVINGREG.DisplacementField,invivo_registered] = imregdemons(tmp_invivo_registered,fixed_GRIN,100,'AccumulatedFieldSmoothing',1.5,'PyramidLevels',1);

cd ([path mouse '\'   num2str(Z) 'um\'])
% 
% 
% clear button
% if exist('anchor_coordiantes_invivo.mat')
%     promptMessage = sprintf('anchors coordinates file exist, to define new anchors?');
%     titleBarCaption = 'Continue?';
%     button = questdlg(promptMessage, titleBarCaption, 'load coordinates', 'define new', 'load coordinates');
%   else
%     button='define new';  
% end
% 
% if strcmpi(button, 'load coordinates')
%     load('anchor_coordiantes_invivo')
% else
%     maxAllowablePoints=3;
%     [x_anchor,y_anchor]=get_cell_coordinates(fixed_GRIN,invivo_registered,maxAllowablePoints);
%     save('anchor_coordiantes_invivo','x_anchor','y_anchor')
% end
%%%%%%% define cells 
% first check if cell file alraedy exist 
ASK=0; % define if open question box or not

clear button
if exist('Cell_coordiantes_invivo.mat') && ASK
    promptMessage = sprintf('cell coordinates file exist, to define new anchors?');
    titleBarCaption = 'Continue?';
    button = questdlg(promptMessage, titleBarCaption, 'load coordinates', 'define new', 'load coordinates');
elseif ~exist('Cell_coordiantes_invivo.mat')
    button='define new';
elseif ~ASK
    button='load coordinates';
end
% load cell coordinates OR open a figure to define 
if strcmpi(button, 'load coordinates')
    load('Cell_coordiantes_invivo')
else
    maxAllowablePoints=40;
    [x_all_cells,y_all_cells]=get_cell_coordinates(fixed_GRIN,invivo_registered,maxAllowablePoints);
    %if (x_all_cells(end)<0) || (y_all_cells(end)<0)
        x_all_cells=x_all_cells(1:end-1);
        y_all_cells=y_all_cells(1:end-1);
    %end
    save('Cell_coordiantes_invivo','x_all_cells','y_all_cells')
end

%% define score threshhold
iteration.alow_iteration=1;
iteration.step=3;
iteration.limit=3;
parameter.radius=radius;
parameter.exponents_array=[0.8 0.8 0.8];
% compare cells using 'ssim' function 
for ci=1:length(x_all_cells)    
    [I2{ci},xc{ci},yc{ci},cell_score(ci),is_cell(ci)]=get_cell_similarity_score(fixed_GRIN,invivo_registered,x_all_cells(ci),y_all_cells(ci),parameter,iteration );
end
%noW plot it
figure;
% plot in-vivo
sh{1}=subplot(1,2,1); imshow(invivo_registered); xlabel('in vivo'); hold on

for ci=1:length(x_all_cells)   
    if is_cell(ci)
        plot(xc{ci},yc{ci},'w--','LineWidth',2); hold on
    end
     th=text(x_all_cells(ci),y_all_cells(ci),num2str(ci),'FontSize',14); hold on
     if cell_score(ci)>=thresh && is_cell(ci)
        th.Color='g';% th.FontSize=18;
     elseif ~is_cell(ci)
        th.Color='k';
     elseif cell_score(ci)<thresh && is_cell(ci)
        th.Color='r';% th.FontSize=18;
     end
end
% plot fixed
sh{2}=subplot(1,2,2); imshow(fixed_GRIN); xlabel('fixed')
hold on
map = my_map('magenta');
sh{2}.Colormap=map;
title([mouse ' '  num2str(Z) 'um'])
for ci=1:length(x_all_cells)   
     if is_cell(ci)
        plot(xc{ci},yc{ci},'w--','LineWidth',2); hold on 
     end
     th=text(x_all_cells(ci),y_all_cells(ci),num2str(ci),'FontSize',14); hold on
     if cell_score(ci)>=thresh && is_cell(ci)
        th.Color='g';% th.FontSize=18;
     elseif ~is_cell(ci)
        th.Color='k';
     elseif cell_score(ci)<thresh && is_cell(ci)
        th.Color='r';% th.FontSize=18;
     end
end
% plot cells
figure
for ci=1:length(x_all_cells)
    subplot(4,ceil(length(x_all_cells)/4),ci)
    montage({I2{ci}{2},I2{ci}{1}}) % 2 is in-vivo. 1 is fixed

    th=title(num2str(ci));
    if cell_score(ci)>=thresh && is_cell(ci)
        th.Color='g';% th.FontSize=18;
        th.FontSize=18;
    elseif ~is_cell(ci)
        th.FontSize=8;
        th.Color='k';
    elseif cell_score(ci)<thresh && is_cell(ci)
        th.Color='r';% th.FontSize=18;
        th.FontSize=18;
    end
end



% get score for shuffeled cells, to check how method is good
r2 = randperm(length(I2),length(I2));
for li=1:length(I2)-1
        ci=r2(li);
        c2i=r2(li+1);
        M_size1=min(size(I2{ci}{1},1),size(I2{c2i}{2},1));
        M_size2=min(size(I2{ci}{1},2),size(I2{c2i}{2},2));
        shuffeled_score(ci)=ssim(I2{ci}{1}(1:M_size1,1:M_size2),I2{c2i}{2}(1:M_size1,1:M_size2),'Exponents',parameter.exponents_array);

end
Last_ind=r2(end);
M_size1=min(size(I2{Last_ind}{1},1),size(I2{1}{2},1));
M_size2=min(size(I2{Last_ind}{1},2),size(I2{1}{2},2));
shuffeled_score(Last_ind)=ssim(I2{Last_ind}{1}(1:M_size1,1:M_size2),I2{1}{2}(1:M_size1,1:M_size2),'Exponents',parameter.exponents_array);



image_score=ssim(invivo_registered,fixed_GRIN,'Exponents',parameter.exponents_array);%luminance, contrast, and structural terms,
peaksnr = psnr(invivo_registered,fixed_GRIN) ; %calculates the peak signal-to-noise ratio for the image
err = immse(invivo_registered,fixed_GRIN);% calculates the mean-squared error (MSE) between the arrays X and Y.
