function Plot_RewardConditioning_processed
%% Anat Kahan, Ryan Cho, Cell Reports 2021
clear all
%mouse='DAT9R'; 
mouse='DAT11RL'; 


Fs=381.47;
fs=382;
before=20;
after=40;
switch mouse
    case 'DAT9R' % this mouse had higher signal 
        filename='D:\DATA_Glab\Behavioral\DAT_water_restricted\DAT9R_061318_processed.mat';
        load(filename)
        x1=Z;
        filename='D:\DATA_Glab\Behavioral\DAT_water_restricted\DAT9R_061418_processed.mat';
        load(filename)
        x2=Z;
        filename='D:\DATA_Glab\Behavioral\DAT_water_restricted\DAT9R_061518_processed.mat';
        load(filename)
        x3=Z;
        
        % exclude data which is mis- aligned
        x1.dF_onset=[x1.dF_onset(1,:)' x1.dF_onset(3:end,:)']';
        x1.lick_onset=[x1.lick_onset(1,:)' x1.lick_onset(3:end,:)']';
        x3.dF_onset=[x3.dF_onset(1:5,:)']';
        x3.lick_onset=[x3.lick_onset(1:5,:)']';
        num_of_cells=283;% based on LiGS, countification using Imaris 9.2.0 'spots'
        distance_factor=5.967312608;% 'distance_factor' based on (exp(-dz)) found on LiGS histology  11RL_5X.xls
    case 'DAT11RL'% this mouse had lower signal 
        filename='D:\DATA_Glab\Behavioral\DAT_water_restricted\DAT11RL_061318_processed.mat';
        load(filename)
        x1=Z;
        filename='D:\DATA_Glab\Behavioral\DAT_water_restricted\DAT11RL_061418_processed.mat';
        load(filename)
        x2=Z;
        filename='D:\DATA_Glab\Behavioral\DAT_water_restricted\DAT11RL_061518_processed.mat';
        load(filename)
        x3=Z;
        
        x1.dF_onset=[x1.dF_onset(1:2,:)' x1.dF_onset(5:end,:)']';
        x2.dF_onset=[x2.dF_onset(1:2,:)' x2.dF_onset(4,:)' x2.dF_onset(6:end,:)']';
        x1.lick_onset=[x1.lick_onset(1:2,:)' x1.lick_onset(5:end,:)']';
       x2.lick_onset=[x2.lick_onset(1:2,:)' x2.lick_onset(4,:)' x2.lick_onset(6:end,:)']';
        num_of_cells=179; % based on LiGS, countification using Imaris 9.2.0 'spots'
        distance_factor=0.942667789; % 'distance_factor' based on (exp(-dz)) found on LiGS histology  11RL_5X.xls
end

c = copper;

all_dF_onset=[x1.dF_onset' x2.dF_onset' x3.dF_onset' ];% combine all trials
all_dF_onset_original=all_dF_onset;
baseline=mad(all_dF_onset(fs*2:7000,:));
all_dF_onset = (all_dF_onset - median(all_dF_onset))./mad(baseline); % normalization using robust z-score    

%all_dF_onset=zscore(all_dF_onset)-mad(all_dF_onset(1:7000,:));
all_lick_onset=[x1.lick_onset' x2.lick_onset' x3.lick_onset' ];
%% SEM calculation
sem_dF=(std(all_dF_onset')/sqrt(size(all_dF_onset,2)))';
sem_dF_org=(std(all_dF_onset_original')/sqrt(size(all_dF_onset_original,2)))';
sem_lick=(std(all_lick_onset')/sqrt(size(all_lick_onset,2)))';

%%
t = -before:1/fs:after;

% %% plot all repeats
% figure
% plot(t,all_dF_onset)
% hold on
% plot(t,all_lick_onset)
% title(mouse)
% xlabel('Time (sec)');
% ylabel('\DeltaF/F (z-score)');
% ylim([-5,23])
num_sec=1+5;
ind_sec=intersect(find(t>1),find(t<num_sec)); % X first seconds of response
disp (['averaged df/f response of ' mouse ' during the first ' num2str(num_sec) ' seconds: ' num2str(mean(mean(all_dF_onset_original(ind_sec,:)))) '+-' num2str(std(mean(all_dF_onset_original(ind_sec,:)))/sqrt(length(mean(all_dF_onset_original(ind_sec,:)))))])
disp (['averaged z-scored df/f response of ' mouse ' during the first ' num2str(num_sec) ' seconds: ' num2str(mean(mean(all_dF_onset(ind_sec,:)))) '+-' num2str(std(mean(all_dF_onset(ind_sec,:)))/sqrt(length(mean(all_dF_onset(ind_sec,:)))))])
disp (['averaged df/f response of ' mouse ' during the first ' num2str(num_sec) ' seconds, norm to distances: ' num2str(mean(mean(all_dF_onset_original(ind_sec,:)))/distance_factor) '+-' num2str(std(mean(all_dF_onset_original(ind_sec,:)))/sqrt(length(mean(all_dF_onset_original(ind_sec,:))))/distance_factor)])

%% plot mean
figure
subplot (1,4,1);
plot(t,mean(all_lick_onset'),'color',[0 0 0])
hold on
plot(t,mean(all_lick_onset')+sem_lick','-','color',[0.5 0.5 0.5],'linewidth',1)
hold on
plot(t,mean(all_lick_onset')-sem_lick','-','color',[0.5 0.5 0.5],'linewidth',1)
ylabel('Licks');
ylim([-2/22*10,8])
xlabel('Time (sec)');
title(mouse)
%yyaxis left plot original df/f
subplot (1,4,2);
plot(t,mean(all_dF_onset_original'),'color',c(end,:) ,'linewidth',5)
hold on
plot(t,mean(all_dF_onset_original')+sem_dF_org','-','color',c(32,:),'linewidth',1)
hold on
plot(t,mean(all_dF_onset_original')-sem_dF_org','-','color',c(32,:),'linewidth',1)
ylabel('\DeltaF/F');
ylim([-5,30])
xlabel('Time (sec)');
title(mouse)
%fill([t; mean(all_dF_onset')-sem_dF'],[t;mean(all_dF_onset')+sem_dF'],[1,0,0],'LineStyle','none')
hold on
%% z-scored
subplot (1,4,3);
plot(t,mean(all_dF_onset'),'color',c(end,:) ,'linewidth',5)
hold on
plot(t,(mean(all_dF_onset')+sem_dF'),'-','color',c(32,:),'linewidth',1)
hold on
plot(t,(mean(all_dF_onset')-sem_dF'),'-','color',c(32,:),'linewidth',1)
ylabel('\DeltaF/F (z-score)');
ylim([-5,30])
xlabel('Time (sec)');
title(mouse)
%% norm to distnace
subplot (1,4,4);
plot(t,mean(all_dF_onset_original')/distance_factor,'color',c(end,:) ,'linewidth',5)
hold on
plot(t,(mean(all_dF_onset_original')+sem_dF_org')/distance_factor,'-','color',c(32,:),'linewidth',1)
hold on
plot(t,(mean(all_dF_onset_original')-sem_dF_org')/distance_factor,'-','color',c(32,:),'linewidth',1)
ylabel('\DeltaF/F (normalized to distance)');
ylim([-1,7])
xlabel('Time (sec)');
title(mouse)

%% plot map
figure
% plot licks
ax2=subplot (1,2,1);
ph=imagesc(all_lick_onset');
ph_xlabel=get(ph.Parent,'xticklabels');
ph_x=get(ph.Parent,'xtick');
%set(ph.Parent,'xticklabels',pixel2time(ph_xlabel,length(all_dF_onset),before,after,fs))
[y_str y_num]=pixel2time(ph_xlabel,length(all_dF_onset),before,after,fs);
set(ph.Parent,'xtick',ph_x)
set(ph.Parent,'xticklabels',y_str)
caxis([0 8])
ylim([0.5 21.5])
title ([mouse ' Licks'])
ylabel('trials');
xlabel('Time (sec)');
colormap(ax2,bone)
colorbar
% plot df
ax1=subplot (1,2,2);
ph=imagesc(all_dF_onset_original');
%text(1000,3,mouse,'FontSize',14)
ph_xlabel=get(ph.Parent,'xticklabels');
ph_x=get(ph.Parent,'xtick');
[y_str y_num]=pixel2time(ph_xlabel,length(all_dF_onset),before,after,fs);
set(ph.Parent,'xtick',ph_x)
set(ph.Parent,'xticklabels',y_str)
caxis([-0 16])
%caxis([-5 23])
ylim([0.5 21.5])
colormap(ax1,copper)
colorbar
title ([mouse ' \DeltaF/F '])
ylabel('trials');
xlabel('Time (sec)');
% % plot df
% ax1=subplot (1,4,3);
% ph=imagesc(all_dF_onset');
% %text(1000,3,mouse,'FontSize',14)
% ph_xlabel=get(ph.Parent,'xticklabels');
% ph_x=get(ph.Parent,'xtick');
% [y_str y_num]=pixel2time(ph_xlabel,length(all_dF_onset),before,after,fs);
% set(ph.Parent,'xtick',ph_x)
% set(ph.Parent,'xticklabels',y_str)
% caxis([-0 35])
% %caxis([-5 23])
% ylim([0.5 21.5])
% colormap(ax1,copper)
% colorbar
% title ([mouse ' \DeltaF/F '])
% ylabel('trials');
% xlabel('Time (sec)');
% % plot df
% ax1=subplot (1,4,4);
% ph=imagesc((all_dF_onset_original./distance_factor)');
% %text(1000,3,mouse,'FontSize',14)
% ph_xlabel=get(ph.Parent,'xticklabels');
% ph_x=get(ph.Parent,'xtick');
% [y_str y_num]=pixel2time(ph_xlabel,length(all_dF_onset),before,after,fs);
% set(ph.Parent,'xtick',ph_x)
% set(ph.Parent,'xticklabels',y_str)
% caxis([-0 4])
% %caxis([-5 23])
% ylim([0.5 21.5])
% colormap(ax1,copper)
% colorbar
% title ([mouse ' \DeltaF/F '])
% ylabel('trials');
% xlabel('Time (sec)');
end

function [y_str y_num] = pixel2time(X,L,before,after,fs)
t = -before:1/fs:after;
dur=after+before+1/fs;
for i=1:length(X)
    y_str{i}=num2str((str2num(X{i})/length(t)*dur*10000-before),'%.0f');
    y_num(i)=(str2num(X{i})/length(t)*dur*10000-before);
    %y{i}=num2str(round((str2num(X{i})*length(t)/L)*dur-before)/3-20);
end

end