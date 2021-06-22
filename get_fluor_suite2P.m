function [norm_all_spks,all_spks,norm_all_dF,all_dF,all_F,cell_ind,cell_score] = get_fluor_suite2P(mouse_info,sess,exp,Fig,Fs)
%% after running suite2P, a Fall.mat file is generated and get fluorscence of 'real' cells 
    % this function plot only the 'cell' , but returns fluorscence of all
%Fig=0;  
clear M cell_ind cell_F norm_cell_F
path='C:\Users\anatk\Documents\Light_sectioning\';
%mouse='WT2170' ;M_num=2170; Z=960; fixed_section=27; fixed_path='WT2170_processed_confocal\WT2170_GC488_ARC_647_10X_zoom2p5_3um_step_10X.tiff_files\';
mouse=mouse_info.mouse ;
%M_num=mouse_info.M_num; Z=mouse_info.Z; fixed_section=mouse_info.fixed_section; fixed_path=mouse_info.fixed_path; 
switch exp
    case 'behavior'
        switch mouse_info.session_style
            case 'combined'
                cd ([path mouse '\all_tiff_to_align_combined\sess1\suite2p\plane0\']) % everything is in one session
            case 'individual'
                switch mouse_info.suite2P_mode
                    case 'sparse'
                        cd ([path mouse '\all_tiff_to_align\sess' num2str(sess) '\suite2p\plane0\'])
                    case 'non_sparse'
                        cd ([path mouse '\all_tiff_to_align_nonSparseMode\sess' num2str(sess) '\suite2p\plane0\'])
                end
        end
    case 'Zstack'
        cd ([path mouse '\tiff_Zstack_to_align\sess' num2str(sess) '\suite2p\plane0\'])
end
%fixedimage=imread([path mouse '\' num2str(Z) 'um\section_' num2str(fixed_section) '_Tissue2GRIN940nm.tif']);

%% load matrix (M) with suite 2P info 
M=load('Fall.mat');

cell_ind=find(M.iscell(:,1)>0);
cell_score=M.iscell(:,2);
%cell_F=M.F(cell_ind,:);
all_F=M.F-0.7*M.Fneu; % correct for background

% switch mouse_info.session_style
%     case 'combined'
%       all_sess_index_on=int64([1:size(all_F,2)/mouse_info.exp_per_sess:size(all_F,2)]);
%       all_sess_index_off=int64([size(all_F,2)/mouse_info.exp_per_sess:size(all_F,2)/mouse_info.exp_per_sess:size(all_F,2)]);
%       all_F=all_F(:,all_sess_index_on(sess):all_sess_index_off(sess)); 
% end

A=size(all_F,2);
for i=1:size(all_F,1) % calcolate dF-F0/F0
    switch mouse_info.session_style
        case 'combined'
            % this is for 2 sessions that their tiff is processed together 
            all_dF(i,:)=[(all_F(i,1:A/2)-median(all_F(i,1:A/4)))/median(all_F(i,1:A/4)) (all_F(i,A/2:end)-median(all_F(i,A/2+1:3*A/4)))/median(all_F(i,A/2+1:3*A/4))];
           
        case 'individual'
            all_dF(i,:)=(all_F(i,:)-median(all_F(i,1:A/4)))/median(all_F(i,1:A/4));
            % all_dF(i,:)=(all_F(i,:)-min(all_F(i,:)))/min(all_F(i,:));
            % all_dF(i,:)=all_F(i,:);
    end
end
all_spks=M.spks;
for i=1:size(M.F,1)
    norm_all_dF(i,:)=all_dF(i,:)./max(all_dF(i,:));
    norm_all_spks(i,:)=M.spks(i,:)./max(M.spks(i,:));
end

norm_cell_dF=norm_all_dF(cell_ind,:);
norm_cell_spks=norm_all_spks(cell_ind,:);
t_array=[0:60/Fs:60*(size(norm_cell_dF,2)-1)/Fs]/3600;% time array in min
%% figure to show traces map
if Fig
figure
subplot(1,3,1)
imagesc(norm_cell_dF)
title([mouse ' ,sess ' num2str(sess) 'norm F map'])

%% plot show traces 
subplot(1,3,2)
for i=1:size(norm_cell_dF,1)
   plot(t_array,norm_cell_dF(i,:)+i)
   hold on
end
title([mouse ' ,sess ' num2str(sess) ' norm F'])

subplot(1,3,3)
for i=1:size(norm_cell_dF,1)
   plot(t_array,norm_cell_spks(i,:)+i)
   hold on
end
title([mouse ' ,sess ' num2str(sess) ' norm spikes'])

end

