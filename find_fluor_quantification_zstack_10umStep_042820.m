function [z,per_Positive]=find_fluor_quantification_zstack_10umStep_042820(mouse)
%% to run this FUNCTION save as tiff in Zen
%% CAN BE USED individually or with 'get_gfap_over_samples'
%
% use this as a helper to count cells. works well for vipGC and vIPGFP.
% but have to be carefull and observe. works well in ~50-60% of images
clear ALL_ch_images ALL_BW_ch_images all_BW_ch_images x y 


if nargin == 0
    %%% set parameters
 %  mouse='VIPGFP11N'
 % mouse='VIPGFP14RL'
  
 % mouse='VIPGC12LL'  
  %mouse='VIPGC25L'  
%  mouse='VIPGC60N'  
 %mouse='VIPGC62L'
  % mouse='VIPGC102L'
 %   mouse='VIPGC103R'
 %mouse='VIPGC106LL'
   mouse='VIPGC107R'
   %mouse='VIPGC108L'
  % mouse='VIPGC110LL'
  % mouse='VIPGC113L'
   % mouse='VIPGC115L'
 % mouse='VIPGC118RR'
   %mouse='VIPGC119L'
   % mouse='VIPGC122R'
    %mouse='VIPGC123L'
  % mouse='VIPGC128R'
  %  mouse='VIPGC129L'
    %mouse='VIPGC162'
    %mouse='VIPGC174RR'
end

show_correlated_figure=0
show3D_fig=0
scn_size=450;% in microns
image_pixel_size=2048;
image_size=850.19;% in microns
scn_pixel_size=image_pixel_size*scn_size/image_size;
count_GFP=0
center='SCN'% choose to count cells above SCN or above fiber 
%center='fiber'
count_GFAP_pixels=1
if exist('ch0')
    move_files=0;
else
    move_files=1;
end

str_to_find='-2048y0';

my_path='Z:\Anat\Data_Glab_from_laptop';
%my_path='D:\DATA_Glab';
%%% tiff files quality is much better here
new_staining=1; % new staining for GFAP/GFP 
if new_staining
    file_name='LGS_SCN_VIP3';ni=1;
else
    file_name='LGS_SCN_VIP';ni=1;
end
[NUMpar,TXTpar,RAWpar]=xlsread([my_path '\' file_name '.xlsx']);
raw_ind=find(strcmp(TXTpar(:,1),mouse));
disp(['mouse raw is ' num2str(raw_ind)])

folder=[TXTpar{raw_ind,7+ni} TXTpar{raw_ind,8+ni}];
step=RAWpar{raw_ind,9+ni}; % image step. can be found in Zen image info
num_channel=RAWpar{raw_ind,10+ni};
fiber_image_ind=RAWpar{raw_ind,14+ni}-1; 
%scn_range=RAWpar{raw_ind,12};
min_X=RAWpar{raw_ind,16+ni};
GFP_ch=RAWpar{raw_ind,15+ni};
scn_image_ind=RAWpar{raw_ind,11+ni};% index of image scn is best expressed
scn_range=[RAWpar{raw_ind,12+ni}:RAWpar{raw_ind,13+ni}];

type='*.tiff';

Factor=0.0;


cd (folder)
if move_files
    for ci=1:num_channel
        dir_ch=ci-1;
        %% read relevant jpegs and move to designated folder and gointo relevnt folder
        move_files_to_folder(folder,dir_ch,type)
    end
end

%% put all the tiff files from each channel into one matrix
for ci=1:num_channel
    clear all_BW_ch_images all_ch_images
    dir_ch=ci-1;
    cd ([folder '/ch' num2str(dir_ch)])
    listoffiles=dir('*.tiff');;
    for i=1:length(listoffiles)
        % finds the index of the image
        tmp_image=imread((listoffiles(i).name));
        this_image=tmp_image(:,:,2);
        %switch mouse
         %   case 'VIPGFP14RL'
          %      k=str2num(listoffiles(i).name(strfind(listoffiles(i).name,'_10X')+6:strfind(listoffiles(i).name,['c' num2str(dir_ch) 'x0'])-1))+1;
           % otherwise
                k=str2num(listoffiles(i).name(strfind(listoffiles(i).name,'b0v0t0z')+7:strfind(listoffiles(i).name,['c' num2str(dir_ch) 'x0'])-1))+1;
        %end
        all_ch_images(:,:,k)=this_image;
    end
    ALL_ch_images{ci}=all_ch_images;
    %  ALL_BW_ch_images{ci}=all_BW_ch_images;
    cd ../
end

%% define SCN rectangle
switch center
    case 'SCN'
      % cd ../
        if exist('scn_coordiantes.mat')
            scn_coor=open('scn_coordiantes.mat');x=scn_coor.x; y=scn_coor.y;
        else
            ci=1;
            figure
            imshow(ALL_ch_images{ci}(:,:,scn_image_ind))
            title(['choose 2 coordinates of middle '  center ' , starting at the bottom'])
            [x,y] = white_ginput(2); % choose 2 cursors, starting bottom
            save('scn_coordiantes','x','y')
        end
    case 'fiber'
         if exist('fiber_1_coordiantes.mat')
            fibercoor=open('fiber_1_coordiantes.mat');
            x=fibercoor.x; y=fibercoor.y;
        else
            ci=1;
            figure
            imshow(ALL_ch_images{ci}(:,:,fiber_image_ind))
            title(['choose 1 coordinates of middle '  center ])
            [x,y] = white_ginput(1); % choose 2 cursors, starting bottom
            save('fiber_1_coordiantes','x','y')
         end       
end

%% count GFP at channel 1 (0 in Zen)
if count_GFP
   [number_cells]= count_cells(x,y,num_channel,ALL_ch_images,scn_pixel_size,step);
end

%% after making sure that GFP counting is working well and xlsx file is updated, it continues to GFAP quantification 
clear per_Positive
if count_GFAP_pixels

    if show_correlated_figure
    cell_count_fig=1
    [scn_z,scn1_cells,scn2_cells]=read_cell_counts(mouse,cell_count_fig);
    end
    clear scn1 scn2 all_scn1_images all_scn2_images all_ch_images scn1_org scn2_org
    for ci=1:num_channel
        % ci=2;
        % all_ch_images=ALL_BW_ch_images{ci};
        all_ch_images_org=ALL_ch_images{ci};
        clear I2
        for i=1:size(all_ch_images_org,3)
            if (y(1)-scn_pixel_size/2)<0; y1=0; else y1=y(1)-scn_pixel_size/2;end
            if (x(1)-scn_pixel_size/2)<0; x1=0; else x1=x(1)-scn_pixel_size/2;end
            scn1_org(:,:,i) = imcrop(all_ch_images_org(:,:,i),[x1 y1 scn_pixel_size scn_pixel_size]);
            
            this_image=scn1_org(:,:,i);
            % creates BW images for each SCN and each plane seperately
            level = graythresh(this_image);
            this_BW = imbinarize(this_image,level);
            scn1(:,:,i)=this_BW;
            clear level this_BW
            if (y(2)-scn_pixel_size/2)<0; y2=0; else y2=y(2)-scn_pixel_size/2;end
            if (x(2)-scn_pixel_size/2)<0; x2=0; else x2=x(2)-scn_pixel_size/2;end
            scn2_org(:,:,i) = imcrop(all_ch_images_org(:,:,i),[x2 y2 scn_pixel_size scn_pixel_size]);
            this_image=scn2_org(:,:,i);
            level = graythresh(this_image);
            this_BW = imbinarize(this_image,level);
            
            scn2(:,:,i)=this_BW;
            clear level this_BW
        end
        if show3D_fig
        figure
        imshow3D(uint8(scn1).*scn1_org)
        figure
        imshow3D(uint8(scn2).*scn2_org)
        end
        all_scn1_images{ci}=scn1;
        all_scn2_images{ci}=scn2;
        all_scn1_images_org{ci}=scn1_org;
        all_scn2_images_org{ci}=scn2_org;
        
        
        for si=1:2 % for 2 scn  % 1 is bottom, 2 is up. head position right always
            clear a
            for i=1:size(all_scn1_images{ci},3)
                clear stats a a_org
                switch si
                    case 1
                        a=all_scn1_images{ci}(:,:,i);
                        a_org=all_scn1_images_org{ci}(:,:,i);
                        
                    case 2
                        a=all_scn2_images{ci}(:,:,i);
                        a_org=all_scn2_images_org{ci}(:,:,i);
                end
                b=uint8(a).*a_org;
                % calculates the sum of fluorscence, normalized to size
                per_Positive(i,ci,si)=sum(sum(b))/(size(b,1)*size(b,2)); % z step/channels/scn side
                z(i)=step*(i-1);
                
            end
        end
        
    end
   z=z'; 
end
 if show_correlated_figure
figure
subplot(2,1,1)
% scn1
[~,ind_max1]=max(scn1_cells(1,:));
scn_z_max(1)=scn_z(ind_max1);

[AX,H1,H2]=plotyy(z,[per_Positive(:,2,1) nan(length(z),1)],scn_z,[nan(length(scn_z),1) scn1_cells(1,:)'],'line','bar')
set(H1,'Color',[0, 0.4470, 0.9410]) % a
set(H2,'FaceColor',[0, 0.5, 0]) % b
%h2(1).FaceColor =[0, 0.5, 0]; h2(2).FaceColor =[0, 0.4470, 0.9410];h2(3).FaceColor =[0.50, 0.50, 0.50]; 
hold on
plot([scn_z_max(1) scn_z_max(1)],[0 max(max(max(per_Positive)))])
%scn2
subplot(2,1,2)
[~,ind_max2]=max(scn2_cells(1,:));
scn_z_max(2)=scn_z(ind_max2);
[AX,H1,H2]=plotyy(z,[per_Positive(:,2,2) nan(length(z),1)],scn_z,[nan(length(scn_z),1) scn2_cells(1,:)'],'line','bar')
set(H1,'Color',[0, 0.4470, 0.9410]) % a
set(H2,'FaceColor',[0, 0.5, 0]) % b
%h2(1).FaceColor =[0, 0.5, 0]; h2(2).FaceColor =[0, 0.4470, 0.9410];h2(3).FaceColor =[0.50, 0.50, 0.50]; 
hold on
plot([scn_z_max(2) scn_z_max(2)],[0 max(max(max(per_Positive)))])

xlabel('z (um)')
ylabel('integrated intensity')
title(mouse)



 end

 %% inorder to return 647 as ch3 - to fit other imaging 
 if   (strcmp(mouse,'VIPGC118RR')|| strcmp(mouse,'VIPGFP12R'))&& ~new_staining% no GFAP staining    
        tmp=[per_Positive(:,1,:) per_Positive(:,1,:).*nan per_Positive(:,2,:)];
        per_Positive=tmp;
 end     

 
 
 
 
