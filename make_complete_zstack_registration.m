function [outputArg1,outputArg2] = make_complete_zstack_registration
% Anat Kahan, Cell Reports 2021
%This function makes a full registration to all zstack, based on manually
%defined points
%
close all

comp=1.10 %1.1 is invivo vs processed fixed; 2.1 is procesessd, GRIN vs tissue side, 1100nm;

%%% read coordinates from xls file
file_name='LGS registration table'; % 3 is with estrus state refered to proestrus
%file_name='LGS registration table home'; % 3 is with estrus state refered to proestrus
path='D:\Data_Glab\Light_sectioning\';
full_path=[path file_name '.xlsx'];
%[NUM,TXT,RAW]=xlsread(full_path,'VS AK');
[NUM,TXT,RAW]=xlsread(full_path);

%%%%%%%%%%
invivo=[];
mouse='WT36R_LGS';num_reg_points=8;this_mouse=strfind(TXT(:,1),mouse(1:end-4));

% finds this mouse index
for ii=1:length(this_mouse); if ~isempty(this_mouse{ii}); this_mouse_ind=ii; end; end
rows=4;
all_Z=-1*NUM(this_mouse_ind+rows:this_mouse_ind+rows+num_reg_points-1,2);
all_s2=NUM(this_mouse_ind+rows:this_mouse_ind+rows+num_reg_points-1,3);
all_s1=NUM(this_mouse_ind+rows:this_mouse_ind+rows+num_reg_points-1,5);
all_s1100=NUM(this_mouse_ind+rows:this_mouse_ind+rows+num_reg_points-1,7);
all_s3=NUM(this_mouse_ind+rows:this_mouse_ind+rows+num_reg_points-1,9);

%%%%%%%%%%%


inds= strfind(mouse,'_');
if length(inds)>1
    M_num=mouse(inds(1)+1:inds(2)-1);
    mouse_title=[mouse(1:inds(1)-1) ' ' M_num ' ' mouse(inds(2)+1:end)];
else
    M_num=mouse(3:inds(1)-1);
    mouse_title=[mouse(1:2) ' ' M_num ' ' mouse(inds(1)+1:end)];
end

switch comp
    case 1.1
        for zi=1:length(all_Z)
            Z=all_Z(zi);
            if strcmp(invivo,'WT36invivo1')
                load(['D:\DATA_Glab\Light_sectioning\' mouse '\fit_files_invivo1\' M_num '_invivo1p1_' num2str(Z) '_FIT_NRS'],'FIT')
            else
                load(['D:\DATA_Glab\Light_sectioning\' mouse '\fit_files\' M_num '_invivo1p1_' num2str(Z) '_FIT_NRS'],'FIT')
            end
            Rs(zi)=FIT.R;
            if isfield(FIT,'BrF');
                BrF1(zi)=FIT.BrF(1); BrF2(zi)=FIT.BrF(2);
            else
                BrF1(zi)=8; BrF2(zi)=8;
            end
            angles(zi)=FIT.this_angle;
            scales(zi)=FIT.this_scale;
            Tinv=invert(FIT.mytform);
            TXs(zi)=Tinv.T(3,1);
            TYs(zi)=Tinv.T(3,2);
        end
    case 2.1
        %   save(['D:\DATA_Glab\Light_sectioning\' mouse '\' M_num '_processed_1100nm_' num2str(Z) '_tform_' method_fixed ],'mytform')
        load(['D:\DATA_Glab\Light_sectioning\' mouse '\' M_num '_processed_1100nm_' num2str(Z) '_FIT_' method_fixed ],'FIT')
    case 2
        for zi=1:length(all_Z)
            Z=all_Z(zi);
            load(['D:\DATA_Glab\Light_sectioning\' mouse '\fit_files\' M_num '_processed_940nm_' num2str(Z) '_FIT_NRS'],'FIT')
            
            Rs(zi)=FIT.R;
            if isfield(FIT,'BrF');
                BrF1(zi)=FIT.BrF(1); BrF2(zi)=FIT.BrF(2);
            else
                BrF1(zi)=8; BrF2(zi)=8;
            end
            angles(zi)=FIT.this_angle;
            scales(zi)=FIT.this_scale;
            Tinv=invert(FIT.mytform);
            TXs(zi)=Tinv.T(3,1);
            TYs(zi)=Tinv.T(3,2);
            
        end
        
end

%% interpolate the new angles and scales for the intermediate images
% also has to define the right section of fixed tissue imaging - doing that by interpolation as
% well
full_Z=[all_Z(1):-15:all_Z(end)];
full_angles=interp1(all_Z(1:end),angles,full_Z);
full_scales=interp1(all_Z(1:end),scales,full_Z);
full_TXs=interp1(all_Z(1:end),TXs,full_Z);
full_TYs=interp1(all_Z(1:end),TYs,full_Z);
full_s2=interp1(all_Z(1:end),all_s2(1:end),full_Z);full_s2=round(full_s2);
full_s1=interp1(all_Z(1:end),all_s1(1:end),full_Z);full_s1=round(full_s1);
full_BrF1s=interp1(all_Z(1:end),BrF1,full_Z);
full_BrF2s=interp1(all_Z(1:end),BrF2,full_Z);

figure
subplot(1,2,1)
plot(full_Z,full_angles,'*b'); hold on
plot(all_Z(1:end),angles,'*r'); hold on
title('Angles')
subplot(1,2,2)
plot(full_Z,full_scales,'*b'); hold on;
plot(all_Z(1:end),scales,'*r'); hold on
title('Scales')
mean_angle=mean(full_angles);

atpath=TXT(this_mouse_ind+2,4);
atpath_2=TXT(this_mouse_ind+2,6);
invivo_path2=TXT(this_mouse_ind+2,2);

all_fixedImage=[];
all_distorted=[];
show_brightness=1;
%for zi=2:5
for zi=1:length(full_Z)
    close all
    switch  comp
        case 1.1
            disp('load ''processed GRIN lens'' data')
            fixedImage=imread([atpath{:} '\section_' num2str(full_s2(zi)) '.tif']);
            % load in vivo data
            disp('load ''in vivo'' data')
            atpath2=[invivo_path2{:} '\s-' num2str(full_Z(zi)) 'nm_v2_nrMC_mean.tif'];
            original_movingImage=imread(atpath2);
        case 2
            atpath=TXT(this_mouse_ind+2,4);
            disp('load ''processed GRIN lens'' data')
            fixedImage=imread([atpath{:} '\section_' num2str(full_s2(zi))  '.tif']);
            % load processed tissue imaging data
            disp('load ''processed Tissue lens'' data')
            movingImagetmp=imread([atpath_2{:} '\section_' num2str(full_s1(zi)) '.tif']);
            original_movingImage=fliplr(movingImagetmp);
            
    end
    
    
    [AmovingImage, AfixedImage]=adjust_brightness(original_movingImage,fixedImage,full_BrF1s(zi),full_BrF2s(zi),show_brightness);

    % after finding angle, scale, tx and ty, I re-make the matrix
    ss=full_scales(zi)*cos(full_angles(zi));
    sc=full_scales(zi)*sin(full_angles(zi));
    new_Tinv.T=[sc, -ss, 0; ss sc, 0; full_TXs(zi), full_TYs(zi), 1];
    Tinv.T=new_Tinv.T;
    recoverd_mytform=invert(Tinv);
    
   
    distorted_angle=imrotate(AmovingImage,-mean_angle);
    distorted_angle_scale=imresize(distorted_angle,1/full_scales(zi));% fit to 1.1
    
    
  
    figure
    subplot (1,4,1), imshow(AmovingImage),title('adj.B fixed Image ')
    subplot (1,4,2), imshow(distorted_angle_scale), title('adj.B moving Image ');set(gca,'xtick',[0 175 350]);set(gca,'xticklabel',{'0','175','350'});
    subplot (1,4,3), imshow(AfixedImage), title('adj.B fixed Image ');set(gca,'xtick',[0 175 350]);set(gca,'xticklabel',{'0','175','350'});%
    subplot (1,4,4), imshowpair(distorted_angle_scale,AfixedImage), title('adj.B moving Image ');set(gca,'xtick',[0 175 350]);set(gca,'xticklabel',{'0','175','350'});
    
    all_fixedImage=[AfixedImage' all_fixedImage];
    IC=imresize(distorted_angle_scale,[512 512]);
    %figure; montage({distorted_angle_scale,IC});
    all_distorted=[IC' all_distorted];
    switch comp
        case 1.1
            %     saveastiff(AfixedImage,['D:\DATA_Glab\Light_sectioning\' mouse '\invivo2fixedGRIN\fixed\fixed_invivo2GRIN_' num2str(full_Z(zi)) '.tif']);
            %      saveastiff(IC,['D:\DATA_Glab\Light_sectioning\' mouse '\invivo2fixedGRIN\moving\moving_invivo2GRIN_'   num2str(full_Z(zi)) '.tif']);
            saveastiff(AfixedImage,['D:\DATA_Glab\Light_sectioning\' mouse '\invivo2fixedGRIN\fixed\mean_angle\fixed_invivo2GRIN_meanAngle_' num2str(full_Z(zi)) '.tif']);
            saveastiff(IC,['D:\DATA_Glab\Light_sectioning\' mouse '\invivo2fixedGRIN\moving\mean_angle\moving_invivo2GRIN_meanAngle_'   num2str(full_Z(zi)) '.tif']);
            
        case 2
            %      saveastiff(AfixedImage,['D:\DATA_Glab\Light_sectioning\' mouse '\fixedGRINtoTissue940\fixed\fixed_invivo2GRIN_' num2str(full_Z(zi)) '.tif']);
            %      saveastiff(IC,['D:\DATA_Glab\Light_sectioning\' mouse '\fixedGRINtoTissue940\moving\moving_invivo2GRIN_'   num2str(full_Z(zi)) '.tif']);
            saveastiff(AfixedImage,['D:\DATA_Glab\Light_sectioning\' mouse '\fixedGRINtoTissue940\fixed\mean_angle\fixed_tissue2GRIN_meanAngle_' num2str(full_Z(zi)) '.tif']);
            saveastiff(IC,['D:\DATA_Glab\Light_sectioning\' mouse '\fixedGRINtoTissue940\moving\mean_angle\moving_tissue2GRIN_meanAngle_'   num2str(full_Z(zi)) '.tif']);
            
    end
    
end

end

