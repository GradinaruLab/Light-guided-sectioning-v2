function MAX_GFP_distance_fiber=find_fluor_exp_above_fiber(mouse)
%% to run this script save Zen images as tiff
%% this function is used to fund fluorscence above fiber

if nargin == 0
    %%% set parameters
    %mouse='VIPGFP14RL'
    %mouse='VIPGC106LL'
    %mouse='VIPGC113L'
    %  mouse='VIPGC115L'
    %mouse='VIPGC162'
    mouse='VIPGC129L';
   %  mouse='VIPGC122R';
    % mouse='VIPGC128R';
    %mouse='VIPGC62L';
end


RI=1; % RIseeBD is 1.48, but is neglactible
str_to_find='-2048y0';
mypath='D:\DATA_Glab';
mypath='C:\Users\anatk\Documents\Data_Glab_home_work';
%%% tiff files quality is much better here
file_name='LGS_SCN_VIP2';ni=1; % second round seeDB - staining for GFP, so looks better. 
%file_name='LGS_SCN_VIP';ni=1;
[NUMpar,TXTpar,RAWpar]=xlsread([mypath '\' file_name '.xlsx']);
raw_ind=find(strcmp(TXTpar(:,1),mouse));

folder=[TXTpar{raw_ind,7+ni} TXTpar{raw_ind,8+ni}];
step=RAWpar{raw_ind,9+ni}; % image step. can be found in Zen image info
num_channel=RAWpar{raw_ind,10+ni};
fiber_image_ind=RAWpar{raw_ind,14+ni}-1; % index of image scn is best expressed
%scn_range=RAWpar{raw_ind,12};
min_X=RAWpar{raw_ind,16+ni};
GFP_ch=RAWpar{raw_ind,15+ni};


type='*.tiff';
% define the channel we are looking at

Factor=0.0;
move_files=1;

%folder='D:\DATA_Glab\Light_sectioning\WT2170\WT2170_processed_confocal\WT2170_GC488_ARC_647_10X_zoom2p5_3um_step_10X.tiff_files';
cd (folder)
if exist('SCN_to_fiber_distance.mat')>0
    load('SCN_to_fiber_distance.mat')
    %MAX_GFP_distance_fiber=SCN_to_fiber_distance;
else
    
    if move_files
        for ci=1:num_channel
            dir_ch=ci-1;
            %% read relevant jpegs and move to designated folder and gointo relevnt folder
            move_files_to_folder(folder,dir_ch,type)
        end
    end
    
    
    dir_ch=0;% for GFP
    % go to relevant folder
    cd(['ch' num2str(dir_ch)])
    listOfFiles2 = dir('*.tiff');
    for fi=1:length(listOfFiles2)
        if contains(listOfFiles2(fi).name ,['z' num2str(fiber_image_ind) 'c' num2str(dir_ch) 'x0'])
            ind=fi;
        end
    end
    
    tmp_image=imread((listOfFiles2(ind).name)); cd ../
    fiber_image=tmp_image(:,:,2);
    
    
    %% show figure and define coordinates of fiber
    
    brightness_factor=10;
    if exist('fiber_4_coordiantes.mat')>0
        load('fiber_4_coordiantes.mat')
    else
        
        BrF=[brightness_factor,brightness_factor];
        [Afiber_image, ~]=adjust_brightness(fiber_image,fiber_image,BrF(1),BrF(2),0);
        figure
        imshow(Afiber_image)
        title([mouse ', choose 4 coordinates, starting 12, clockwise'])
        [x,y] = white_ginput(4); % choose 4 cursors, clock wise, starting 12
        save('fiber_4_coordiantes','x','y')
    end
    %% using those input, define center and delta
    Dia=0.5*(abs(x(4)-x(2))+abs(y(3)-y(1)));% calculate diameter
    I2 = imcrop(fiber_image,[min(x)-Dia*Factor/2 min(y)-Dia*Factor/2 Dia*(1+Factor) Dia*(1+Factor)]);
    
    figure
    subplot(1,2,1)
    imshow(fiber_image*brightness_factor)
    subplot(1,2,2)
    imshow(I2*brightness_factor)
    title('check if cropping is good')
    
    %% calculates fluorscence intensity in each channel.
    % use show_image to make sure it is in the desired area
    clear fluor_int
    show_image=0;
    for ci=1:num_channel
        dir_ch=ci-1;
        clear I3 this_image tmp_image BW BWI I3_clean
        cd(['ch' num2str(dir_ch)])
        listOfFiles2 = dir('*.tiff');
        %% find background value to substract
        %     tmp_image=imread((listOfFiles2(ind).name));
        %     I4 = imcrop(tmp_image,[min(x)-Dia*Factor/2 min(y)-Dia*Factor/2 Dia*(1+Factor) Dia*(1+Factor)]);
        %     [J, rect] = imcrop(I4);
        %     J=J(:,:,2);
        for fi=1:length(listOfFiles2)
            tmp_image=imread((listOfFiles2(fi).name));
            this_image=tmp_image(:,:,2);
            %         if fi==1 && ci==1
            %             figure;[J, rect] = imcrop(this_image);
            %         else
            %             J=imcrop(this_image,rect);
            %         end
            I3 = imcrop(this_image,[min(x)-Dia*Factor/2 min(y)-Dia*Factor/2 Dia*(1+Factor) Dia*(1+Factor)]);
            level = graythresh(I3);
            BW = imbinarize(I3,level);
            %         level = graythresh(this_image);
            %         BW = imbinarize(this_image,level);
            %         I3 = imcrop(this_image,[min(x)-Dia*Factor/2 min(y)-Dia*Factor/2 Dia*(1+Factor) Dia*(1+Factor)]);
            BWI=uint8(BW);
            if show_image
                figure; subplot(2,2,1); imshow(I3);
                subplot(2,2,2); imshowpair(I3,BW,'montage')
                subplot(2,2,4);imshowpair(I3,BWI.*I3,'montage'); title(['channel ' num2str(ci) ' plane ' num2str(fi)])
            end
            
            I3_clean=BWI.*I3;
            % find the right index (list is not sorted)
            str1=listOfFiles2(fi).name(strfind(listOfFiles2(fi).name,str_to_find)-8:end);
            if fiber_image_ind>1
                %k=str2num(str1(strfind(str1,'z')+1:strfind(str1,'c')-1))+1-fiber_image_ind;
                k=str2num(str1(strfind(str1,'z')+1:strfind(str1,'c')-1))+1;
            else
                k=str2num(str1(strfind(str1,'z')+1:strfind(str1,'c')-1))+1;
                % if k==0; k=k+1; end
            end
            % put into matrix
            fluor_int(ci,k)=sum(sum(I3_clean));
        end
        cd ../
    end
    % create x and normalize fluor_int
    clear norm_fluor_int X
    for ci=1:num_channel
        cd(['ch' num2str(dir_ch)])
        listOfFiles2 = dir('*.tiff');
        dir_ch=ci-1;
        X(ci,:)=[0:step:(length(listOfFiles2)-1)*step];
        norm_fluor_int(ci,:)=fluor_int(ci,:)./max(fluor_int(ci,:));
        cd ../
    end
    %X=mean(X);
    
    
    figure; subplot(2,1,1); h=bar(X(1,:)',norm_fluor_int(1,:)');
    h(1).FaceColor =[0, 0.5, 0];% h(2).FaceColor =[0, 0.4470, 0.9410]; %h(3).FaceColor =[0.50, 0.50, 0.50];
    title(mouse); xlabel('depth (um)'); ylabel ('Normalized dF (a.u.)')
    %xlim([min_X max(mean(X))])
    subplot(2,1,2); h=bar(X(1,:)',fluor_int(1,:)');
    h(1).FaceColor =[0, 0.5, 0]; %h(2).FaceColor =[0, 0.4470, 0.9410]; %h(3).FaceColor =[0.50, 0.50, 0.50];
    xlabel('depth (um)'); ylabel ('dF (a.u.)')
    % xlim([min_X max(mean(X))])
    title([mouse ' ,choose GFP max'])
    [x_max,y_max] = ginput(1); % choose 4 cursors, clock wise, starting 12
    
    1
    
    if size(X,1)>2
        figure; subplot(2,1,1); h=bar(X',norm_fluor_int');
        h(1).FaceColor =[0, 0.5, 0]; h(2).FaceColor =[0, 0.4470, 0.9410]; h(3).FaceColor =[0.5, 0.5, 0.5];
        title(mouse); xlabel('depth (um)'); ylabel ('Normalized dF (a.u.)')
        xlim([min_X max(mean(X))])
        subplot(2,1,2); h=bar(X',fluor_int');
        h(1).FaceColor =[0, 0.5, 0]; h(2).FaceColor =[0, 0.4470, 0.9410]; h(3).FaceColor =[0.5, 0.5, 0.5];
        xlabel('depth (um)'); ylabel ('dF (a.u.)')
        xlim([min_X max(mean(X))])
    end
    
    %[val,ind]=max(fluor_int(GFP_ch,:));
    MAX_GFP_distance_fiber=RI*abs(fiber_image_ind*step-round(x_max,2,'significant'));
    save('SCN_to_fiber_distance','MAX_GFP_distance_fiber')
end