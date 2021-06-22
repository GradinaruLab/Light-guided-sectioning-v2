function [distance_fiber_scn_center] = find_distance_fiber_to_target(mouse)

%% this function is used to find distance of fiber from target, based on coordinates found using 
%'cell_quantification_zstack_assitance' for fiber  coordinates
% and 'find_fluor_quantification_zstack_10umStep_042820' for target
% coordinates


if nargin == 0
    %%% set parameters
    %mouse='VIPGFP14RL'
    mouse='VIPGC106LL'
    %mouse='VIPGC113L'
    %  mouse='VIPGC115L'
    %mouse='VIPGC162'
    %mouse='VIPGC129L';
    % mouse='VIPGC128R';
    %mouse='VIPGC62L';
end


mypath='D:\DATA_Glab';
mypath='C:\Users\anatk\Documents\Data_Glab_home_work';
%%% tiff files quality is much better here
%file_name='LGS_SCN_VIP2';ni=1;;
file_name='LGS_SCN_VIP';ni=1;
[NUMpar,TXTpar,RAWpar]=xlsread([mypath '\' file_name '.xlsx']);
raw_ind=find(strcmp(TXTpar(:,1),mouse));

folder=[TXTpar{raw_ind,7+ni} TXTpar{raw_ind,8+ni}];
step=RAWpar{raw_ind,9+ni}; % image step. can be found in Zen image info
num_channel=RAWpar{raw_ind,10+ni};
fiber_image_ind=RAWpar{raw_ind,14+ni}-1; % index of image scn is best expressed
%scn_range=RAWpar{raw_ind,12};
min_X=RAWpar{raw_ind,16+ni};
GFP_ch=RAWpar{raw_ind,15+ni};
side=RAWpar{raw_ind,5+ni};
switch side; case 'R'; side_ind=2; case 'L'; side_ind=1; end


% load fiber coordinates
cd(folder)
clear x y
load('fiber_4_coordiantes.mat')
Xfiber=abs(x(4)-x(2));
Yfiber=abs(y(3)-y(1));
% load target coordinates

clear x y
load('scn_coordiantes2.mat') % first index is bottom, which is L
Xtarget=x(side_ind);
Ytarget=y(side_ind);

load('SCN_to_fiber_distance')

pixel_to_um=0.42;

distance_fiber_scn_center=sqrt((pixel_to_um*(Xfiber-Xtarget))^2+(pixel_to_um*(Yfiber-Ytarget))^2+MAX_GFP_distance_fiber^2);
