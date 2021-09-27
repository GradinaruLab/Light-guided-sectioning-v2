function run_Zstack_get_fluor_suite2P
% Anat Kahan, Cell Reports 2021
%% %% this function is used after get_data_suite2P, and inspection of cells that were also stained using LGS
% This function is used to compare rates/amplitudes and width of identified
% cells, 
clear norm_all_F all_F cell_ind F exp
exp_ID='36R';
switch exp_ID
      case '36R'
        mouse_info.mouse='WT36R_LGS'; exp='Zstack' ;mouse_info.M_num=36;
      %  mouse_info.Z=[];% 
        all_Z=[395,485,755];
        sess_stained_cells={[4,14,19,20],[7,14,21,37],16};
        Fs=1;% sampling frame rate
end
this_path='C:\Users\anatk\Documents\Light_sectioning\';

per=2;

for si=1:length(all_Z); sess_markers{si}=['o'];end
% cell indexes, by eye. added one due to python indexing to matlab indexing
sess_fig=1;
for si=1:length(all_Z)
    % this function plot only the 'cell' , but returns fluorscence of all
    [norm_all_spks{si}, all_spks{si}, norm_all_F{si},all_F{si},cell_ind{si},cell_score{si}] = get_fluor_suite2P(mouse_info,all_Z(si),exp,sess_fig,Fs);
    arrays_length(si)=size(norm_all_F{si},2);
end
% rule out the stained positive cells
for si=1:length(all_Z)
    cell_ind{si}=setdiff(cell_ind{si},sess_stained_cells{si});
end
%% z-score  and set all same length
clear all_zscored_F all_zscored_spks
for si=1:size(norm_all_F,2)
    clear this_zscore_F
    norm_all_spks{si}=norm_all_spks{si}(:,1:min(arrays_length));% set norm all spks to be the same size
    all_spks{si}=all_spks{si}(:,1:min(arrays_length));% set all spks to be the same size
    norm_all_F{si}=norm_all_F{si}(:,1:min(arrays_length));% set all F to be the same size
    all_F{si}=all_F{si}(:,1:min(arrays_length));% set all F to be the same size
    this_F=all_F{si};
    this_spks=all_spks{si};
    for ci=1:size(this_F,1)
        clear X Y this_zscore_spks this_zscore_F
        X=zscore(this_F(ci,:));
        Y = prctile(X,per);
        this_zscore_F(ci,:)=X-Y;
        
        Xs=zscore(this_spks(ci,:));
        Ys=prctile(Xs,per);
        this_zscore_spks(ci,:)=Xs-Ys;
    end
    all_zscored_F{si}=this_zscore_F;
    all_zscored_spks{si}=this_zscore_spks;
end

%Get the norm fluorescence of the relevant cells
% first for all stained cells:
k=0;
for si=1:length(sess_stained_cells)
    for ci=1:length(sess_stained_cells{si})
        k=k+1;
        Stained_spks_zscored(k,:)=all_zscored_spks{si}(sess_stained_cells{si}(ci),:);
        Stained_F_zscored(k,:)=all_zscored_F{si}(sess_stained_cells{si}(ci),:);
        Stained_spks_norm(k,:)=norm_all_spks{si}(sess_stained_cells{si}(ci),:);
        Stained_F_norm(k,:)=norm_all_F{si}(sess_stained_cells{si}(ci),:);
    end
end
% next for all the non-stained cells
All_spks_zscored=[];
All_Fzscored=[];
All_spks_norm=[];
All_F_norm=[];
for si=1:length(all_zscored_spks)
    All_spks_zscored=cat(1,All_spks_zscored,all_zscored_spks{si}(cell_ind{si},:));
    All_Fzscored=cat(1,All_Fzscored,all_zscored_F{si}(cell_ind{si},:));
    All_spks_norm=cat(1,All_spks_norm,norm_all_spks{si}(cell_ind{si},:));
    All_F_norm=cat(1,All_F_norm,norm_all_F{si}(cell_ind{si},:));   
end
% plot the non stained 
figure
subplot(1,4,1); imagesc(All_spks_zscored)
subplot(1,4,2); imagesc(All_Fzscored)
subplot(1,4,3); imagesc(All_spks_norm)
subplot(1,4,4); imagesc(All_F_norm)
title([mouse_info.mouse ' non stained cells'])
% plot the stained
figure
subplot(1,4,1); imagesc(Stained_spks_zscored)
subplot(1,4,2); imagesc(Stained_F_zscored)
subplot(1,4,3); imagesc(Stained_spks_norm)
subplot(1,4,4); imagesc(Stained_F_norm)
title([mouse_info.mouse ' stained cells'])

% plot just F_norm, stained vs non-stained
figure
ax(1)=subplot(86,1,1:77); imagesc(All_F_norm)
colormap(ax(1),'gray')
colorbar
ylabel('Cells')
title([mouse_info.mouse ' stained vs non-stained cells'])
ax(2)=subplot(86,1,78:86); imagesc(Stained_F_norm)
[map] = my_map('Magenta');
colormap(ax(2),map)
colorbar
xlabel('Time (sec)')


t_array=[0:60/Fs:60*(min(arrays_length)-1)/Fs];% time array in minutes
% plot the stained cells 
figure
for i=3:size(Stained_spks_zscored,1)
subplot(1,4,1); plot(t_array,Stained_spks_zscored(i,:)+i); hold on;
subplot(1,4,2); plot(t_array,Stained_F_zscored(i,:)+i); hold on;
subplot(1,4,3); plot(t_array,Stained_spks_norm(i,:)+i); hold on;
subplot(1,4,4); plot(t_array,Stained_F_norm(i,:)+i); hold on;
end
title([mouse_info.mouse ' stained cells'])



peak_thresh=0.05;
% find all cells peaks
for ci=1:size(All_spks_zscored,1)  
       %  figure
        %findpeaks(All_F_norm(ci,:),t_array,'Annotate','extents','MinPeakProminence',peak_thresh);
        [allpks{ci},alllocs{ci},allw{ci},allp{ci}]=findpeaks(All_F_norm(ci,:),t_array,'Annotate','extents','MinPeakProminence',peak_thresh);
      %  is_burst(alllocs{si}{ci})       
        rates(ci)=length(allpks{ci})/((min(arrays_length)-1)/Fs);% num peaks per minute 
        peak_proms_mean(ci)=mean(allp{ci});
        peak_width_mean(ci)=mean(allw{ci});       
end
% find stained cells peaks
for ci=1:size(Stained_F_norm,1)  
       %  figure
        %findpeaks(All_F_norm(ci,:),t_array,'Annotate','extents','MinPeakProminence',peak_thresh);
        [Spks{ci},Slocs{ci},Sw{ci},Sp{ci}]=findpeaks(Stained_F_norm(ci,:),t_array,'Annotate','extents','MinPeakProminence',peak_thresh);
      %  is_burst(alllocs{si}{ci})       
        S_rates(ci)=length(Spks{ci})/((min(arrays_length)-1)/Fs);% num peaks per minute 
        S_peak_proms_mean(ci)=mean(Sp{ci});
        S_peak_width_mean(ci)=mean(Sw{ci});       
end

%%%% compare non stained to stained cells
trait_names={'rates', 'amplitude','width'};

for ti=1:length(trait_names)
    trait=trait_names{ti};
    switch trait
        case 'rates'
    specific_cells_trait{1}=rates;  specific_cells_trait{2}=S_rates;
        case 'amplitude'
             specific_cells_trait{1}=peak_proms_mean;  specific_cells_trait{2}=S_peak_proms_mean;
        case 'width'
             specific_cells_trait{1}=peak_width_mean;  specific_cells_trait{2}=S_peak_width_mean;
    end
    sess_names_trait{1}='all cells'; sess_names_trait{2}='stained cells';
     [hsc(ti),psc(ti)]=kstest2(specific_cells_trait{1},specific_cells_trait{2});

    figure
    plotSpread(specific_cells_trait, ...
        'xNames', sess_names_trait, ...
        'distributionMarkers', 'o','distributionColors','r');
    title(trait_names{ti})
end
    
   