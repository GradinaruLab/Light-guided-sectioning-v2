%%cell counting and integrated fluor inyensity quantification

clear ALL_scn_sides z per_Positive
GFAP_OVX=1;
GFAP_count=0;
new_quantifiction=0;
%%% set data set: 
% ALL_mouse={'VIPGFP11N','VIPGFP12R','VIPGFP14RL','VIPGC106LL', 'VIPGC107R', 'VIPGC108L', 'VIPGC110LL', 'VIPGC113L','VIPGC115L','VIPGC118RR', 'VIPGC119L', 'VIPGC122R', 'VIPGC123L', 'VIPGC128R', 'VIPGC129L'};
% % scn_choosen_ind are determined by observation using Zen 
% ALL_mouse_scn_choosen_ind=[22 23 24;21 26 31;10 14 17;14 15 16;8 11 14;14 21 31;23 24 26;13 16 19;20 24 29;15 16 18;8 11 15; 12 16 17;23 26 29;16 17 19;11 12 13];
% ALL_genders={'F','F','F','M','F_OVX','F','F','M','F_OVX','F','F_OVX','F_OVX','F_OVX','F_OVX','NA'};
% ALL_types={'SCN','NA','SCN','SCN','F_OVX','SCN','SCN','SCN_damage','F_OVX','NA','F_OVX','F_OVX','F_OVX','F_OVX','NA'};
% ALL_staining={'Esr1','Esr1','PR','PR','Esr1','PR','PR','PR','Esr1','Esr1','Esr1','PR','PR','PR','Esr1'};
% %recorded_SCN={'R','R','R','R','R','L','R','R','L','R','R','R','R','R'};% 110L 'R by image

ALL_mouse={'VIPGFP11N','VIPGC103R','VIPGC106LL', 'VIPGC107R', 'VIPGC108L', 'VIPGC110LL', 'VIPGC113L','VIPGC115L','VIPGC118RR', 'VIPGC119L', 'VIPGC122R', 'VIPGC123L', 'VIPGC128R', 'VIPGC129L','VIPGC174RR'};

% scn_choosen_ind are determined by observation using Zen 
ALL_mouse_scn_choosen_ind=[22 23 24;19 27 34;14 15 16;8 14 20;14 21 31;23 24 26;13 16 19;20 24 29;15 16 18;8 11 15; 12 16 17;23 26 29;16 17 19;11 12 13;16 19 23];
ALL_genders={'F','M','M','F_OVX','F','F','M','F_OVX','F','F_OVX','F_OVX','F_OVX','F_OVX','NA','F'};
ALL_types={'SCN','SCN','SCN','F_OVX','SCN','SCN','SCN_damage','F_OVX','NA','F_OVX','F_OVX','F_OVX','F_OVX','NA','SCN'};
ALL_staining={'Esr1','Esr1','PR','Esr1','PR','PR','PR','Esr1','Esr1','Esr1','PR','PR','PR','Esr1','PR'};
%recorded_SCN={'R','R','R','R','R','L','R','R','L','R','R','R','R','R'};% 110L 'R by image
recorded_SCN={'nan','R','R','nan','R','R','nan','nan','nan','nan','nan','nan','nan','R','R'};% L is bottom. determined by image 


if GFAP_OVX
    ALL_mouse={'VIPGFP11N', 'VIPGC107R', 'VIPGC108L', 'VIPGC110LL','VIPGC115L','VIPGC118RR', 'VIPGC119L', 'VIPGC122R', 'VIPGC123L', 'VIPGC128R','VIPGC174RR'};
    % scn_choosen_ind are determined by observation using Zen
    ALL_mouse_scn_choosen_ind=[22 23 24;8 11 14;14 21 31;23 24 26;20 24 29;24 30 35;8 11 15; 12 16 17;23 26 29;16 17 19;16 19 23];
    ALL_genders={'F','F_OVX','F','F','F_OVX','F','F_OVX','F_OVX','F_OVX','F_OVX','F'};
    ALL_types={'SCN','F_OVX','SCN','SCN','F_OVX','NA','F_OVX','F_OVX','F_OVX','F_OVX','SCN'};
    ALL_staining={'Esr1','Esr1','PR','PR','Esr1','Esr1','Esr1','PR','PR','PR','PR'};
    %recorded_SCN={'R','R','R','R','R','L','R','R','L','R','R','R','R','R'};% 110L 'R by image
    recorded_SCN={'nan','nan','R','R','nan','nan','nan','nan','nan','nan','R'};% L is bottom. determined by image
end



if GFAP_count
    ALL_mouse={'VIPGC103R','VIPGC106LL', 'VIPGC110LL', 'VIPGC129L','VIPGC174RR'};
    % scn_choosen_ind are determined by observation using Zen
    ALL_mouse_scn_choosen_ind=[19 27 34;14 15 16;14 21 31;11 12 13;16 19 23];
    ALL_genders={'M','M','F','F','F'};
    ALL_types={'SCN','SCN','SCN','SCN','SCN'};
    ALL_staining={'Esr1','PR','PR','Esr1','PR'};
    %recorded_SCN={'R','R','R','R','R','L','R','R','L','R','R','R','R','R'};% 110L 'R by image
    recorded_SCN={'R','R','R','R','R'};% L is bottom. determined by image
end
% my_path='C:\Users\anatk\Documents\Data_Glab_home_work';
% %my_path='D:\DATA_Glab';
% %%% tiff files quality is much better here
% file_name='LGS_SCN_VIP3';ni=1;
% [NUMpar,TXTpar,RAWpar]=xlsread([my_path '\' file_name '.xlsx']);

% this indexing is to compare GFAP recorded side to non recorded in non-OVX
% females/males


ALL_SCN_int=[];

%% define types of gender/ovx with numbers
for i=1:length(ALL_genders)
    if strcmp(ALL_staining{i},'Esr1');       ALL_staining_num(i)=1;
    elseif strcmp(ALL_staining{i},'Esr2');   ALL_staining_num(i)=2;
    elseif strcmp(ALL_staining{i},'PR')      ALL_staining_num(i)=3;
    end
end
%% define types of gender/ovx with numbers
for i=1:length(ALL_genders)
    if strcmp(ALL_genders{i},'F');         ALL_genders_num(i)=1;
    elseif strcmp(ALL_genders{i},'F_OVX'); ALL_genders_num(i)=2;
    elseif strcmp(ALL_genders{i},'M')      ALL_genders_num(i)=3;
    elseif strcmp(ALL_genders{i},'NA')     ALL_genders_num(i)=4;
    end
end
%% define types with numbers
% for i=1:length(ALL_types)
%     if strcmp(ALL_types{i},'SCN');            ALL_types_num(i)=1;
%     elseif strcmp(ALL_types{i},'F_OVX');      ALL_types_num(i)=2;
%     elseif strcmp(ALL_types{i},'SCN_damage'); ALL_types_num(i)=3;
%     elseif strcmp(ALL_types{i},'NA');         ALL_types_num(i)=4;   
%     end
% end
%% define scn side by numbers. 1 is bottom, Left. 2 upper, whish is Right  
for i=1:length(recorded_SCN)
    if strcmp(recorded_SCN{i},'L');         ALL_scn_sides(i)=1; scn_not_recorded(i)=2;
    elseif strcmp(recorded_SCN{i},'R');     ALL_scn_sides(i)=2; scn_not_recorded(i)=1;
    elseif strcmp(recorded_SCN{i},'nan');     ALL_scn_sides(i)=0; scn_not_recorded(i)=0;
    end
end

ALL_genders_num;


if new_quantifiction
    clear z per_Positive 
    for i=1:length(ALL_mouse)
        [z{i},per_Positive{i}]=find_fluor_quantification_zstack_10umStep_042820(ALL_mouse{i});% per positive index: step/channels/scn side
    end
    %Z:\Anat\Data_Glab_from_laptop\Histology\Rep_SCN
    if GFAP_count
        save('Z:\Anat\Data_Glab_from_laptop\Histology\Rep_SCN\SCN_LGS_fluor_quantification3.mat','z','per_Positive');
    elseif GFAP_OVX
        save('Z:\Anat\Data_Glab_from_laptop\Histology\Rep_SCN\SCN_LGS_fluor_quantification4.mat','z','per_Positive');  
    else
      save('Z:\Anat\Data_Glab_from_laptop\Histology\Rep_SCN\SCN_LGS_fluor_quantification2.mat','z','per_Positive');  
    end
 else
     if GFAP_count
         load ('Z:\Anat\Data_Glab_from_laptop\Histology\Rep_SCN\SCN_LGS_fluor_quantification3.mat');
     elseif GFAP_OVX
         load('Z:\Anat\Data_Glab_from_laptop\Histology\Rep_SCN\SCN_LGS_fluor_quantification4.mat');
         
     else
         load('Z:\Anat\Data_Glab_from_laptop\Histology\Rep_SCN\SCN_LGS_fluor_quantification2.mat');
     end
end
%%%%%%%%%
%% arrange integrated intensities data for figure
% reduce data demaintions by averge over 3 slices of imaging
k=0;
m=0;
clear  per_positive_mean per_positive_recorded per_positive_NOT_recorded genders staining those_mice
for i=1:length(z) 
    if ~isempty(z{i})
        k=k+1;
        genders(k)=ALL_genders_num(i);
        staining(k)=ALL_staining_num(i);
        those_mice{k}=ALL_mouse{i};
       % this_z=z{i}(ALL_mouse_scn_choosen_ind(i,:));
        this_mouse_per_positive=mean(per_Positive{i}(ALL_mouse_scn_choosen_ind(i,:),:,:),1);
        %lims=[ALL_mouse_scn_choosen_ind(i,1):ALL_mouse_scn_choosen_ind(i,3)];
        %this_mouse_per_positive=median(per_Positive{i}(lims,:,:),1);
        
        % first index is mouse, second is channel, third is scn side 
        % I choose mean value for 
        per_positive_mean(k,:)=mean(this_mouse_per_positive(1,:,:),3);
        if ALL_scn_sides(i)>0
            m=m+1;
            per_positive_recorded(m,:)=this_mouse_per_positive(1,:,ALL_scn_sides(i));
            per_positive_NOT_recorded(m,:)=this_mouse_per_positive(1,:,scn_not_recorded(i));
        end
    end
end



Esr1_F=intersect(find(staining==1),find(genders==1));%Esr1, female
Esr1_FOVX=intersect(find(staining==1),find(genders==2));%Esr1, female OVX
those_mice(Esr1_FOVX)
PR_F=intersect(find(staining==3),find(genders==1));%PR, female
PR_FOVX=intersect(find(staining==3),find(genders==2));%PR, female OVX
those_mice(PR_FOVX)



combined=[per_positive_recorded; per_positive_NOT_recorded];
% norm_per_positive_recorded=[per_positive_recorded(:,1,:)./max(combined(:,1,:)) per_positive_recorded(:,2,:)./max(combined(:,2,:)) per_positive_recorded(:,3,:)./max(combined(:,3,:))];
% norm_per_positive_NOT_recorded=[per_positive_NOT_recorded(:,1,:)./max(combined(:,1,:)) per_positive_NOT_recorded(:,2,:)./max(combined(:,2,:)) per_positive_NOT_recorded(:,3,:)./max(combined(:,3,:))];

% relative norm
norm_per_positive_recorded=[per_positive_recorded(:,1,:)./(per_positive_recorded(:,1,:)+per_positive_NOT_recorded(:,1,:)) per_positive_recorded(:,2,:)./(per_positive_recorded(:,2,:)+per_positive_NOT_recorded(:,2,:))  per_positive_recorded(:,3,:)./(per_positive_recorded(:,3,:)+per_positive_NOT_recorded(:,3,:)) ];
norm_per_positive_NOT_recorded=[per_positive_NOT_recorded(:,1,:)./(per_positive_recorded(:,1,:)+per_positive_NOT_recorded(:,1,:)) per_positive_NOT_recorded(:,2,:)./(per_positive_recorded(:,2,:)+per_positive_NOT_recorded(:,2,:))  per_positive_NOT_recorded(:,3,:)./(per_positive_recorded(:,3,:)+per_positive_NOT_recorded(:,3,:)) ];

% Z-score normalization 
% norm_per_positive_recorded=[(per_positive_recorded(:,1,:)-sum(combined(:,1,:)))./median(combined(:,1,:)) (per_positive_recorded(:,2,:)-median(combined(:,2,:)))./median(combined(:,2,:)) (per_positive_recorded(:,3,:)-median(combined(:,3,:)))./median(combined(:,3,:))];
% norm_per_positive_NOT_recorded=[(per_positive_NOT_recorded(:,1,:)-median(combined(:,1,:)))./median(combined(:,1,:)) (per_positive_NOT_recorded(:,2,:)-median(combined(:,2,:)))./median(combined(:,2,:)) (per_positive_NOT_recorded(:,3,:)-median(combined(:,3,:)))./median(combined(:,3,:))];

if GFAP_count
    %% GFAP/GFP in implant side vs opposite
    
    figure
    subplot(1,2,1)
    GFP_ind=1;
    GFP_norm=[norm_per_positive_recorded(:,GFP_ind) norm_per_positive_NOT_recorded(:,GFP_ind)];
    ph=bar ([1 2],[mean(GFP_norm(:,1)) mean(GFP_norm(:,2))]);
    ph.FaceColor=[0.5 0.5 0.5];hold on % OVX
    hold on
    plot(GFP_norm','-ok')
    ylim([0 1.1])
    xticks([1 2])
    ylabel('mean integrated fluorescence')
    xticklabels({'Implant side GFP' 'Opposite side GFP'})
    
    subplot(1,2,2)
    GFAP_ind=2;
    GFAP_norm=[norm_per_positive_recorded(:,GFAP_ind) norm_per_positive_NOT_recorded(:,GFAP_ind)];
    ph=bar ([1 2],[mean(GFAP_norm(:,1)) mean(GFAP_norm(:,2))]);
    ph.FaceColor=[0.5 0.5 0.5];hold on % OVX
    hold on
    plot(GFAP_norm','-ok')
    xticks([1 2])
    ylabel('mean integrated fluorescence')
    xticklabels({'Implant side GFAP' 'Opposite side GFAP'})
    ylim([0 1.1])
    
    
    
    % statistics
    h3(1)=lillietest(GFP_norm(:,1));
    h3(2)=lillietest(GFP_norm(:,2));
    if sum(h3)==2
        [HscnG,PscnG]=ttest2(GFP_norm(:,1),GFP_norm(:,2));
    else
        %[PscnG,HscnG]=ranksum(GFP_norm(:,1),GFP_norm(:,2));
        [PscnG,tbl,stats]=kruskalwallis([GFP_norm(:,1),GFP_norm(:,2)]);
    end
    h4(1)=lillietest(GFAP_norm(:,1));
    h4(2)=lillietest(GFAP_norm(:,2));
    
    
    if sum(h4)==2
        [HscnGF,PscnGF]=ttest2(GFAP_norm(:,1),GFAP_norm(:,2));
    else
        [PscnGF,HscnGF]=ranksum(GFAP_norm(:,1),GFAP_norm(:,2));
        [PscnGF,tblGF,statsGF]=kruskalwallis([GFAP_norm(:,1),GFAP_norm(:,2)]);
    end
    
    
    disp(ALL_mouse(ALL_scn_sides>0))
    disp(ALL_genders(ALL_scn_sides>0))
    
    
    %% an alternative figure, based on conversation with Justin Bois, May 2021
    
    figure
    subplot(1,2,1)
    plot(GFAP_norm(:,1)',ones(1,length(GFAP_norm(:,1))),'ok'); hold on
    plot(GFP_norm(:,1)',2*ones(1,length(GFP_norm(:,1))),'ok'); hold on
    plot([0.5 0.5],[0 3],'-'); hold on
    ylim([0 3])
    xlim([0 1])
    xlabel('Relative F, implanted')
    
    
    subplot(1,2,2)
    F=GFAP_norm(:,1)./GFP_norm(:,1);
    plot(F,ones(1,length(F)),'ok'); hold on
    plot([1 1],[0 2],'-'); hold on
     xlim([0.7 2])
     ylim([0 2])
     
     
       h5(1)=lillietest(GFAP_norm(:,1));
    h5(2)=lillietest(GFP_norm(:,1));
    
    
    if sum(h5)==2
        [HscnI,PscnI]=ttest2(GFAP_norm(:,1),GFP_norm(:,1));
    else
        [PscnI,HscnI]=ranksum(GFAP_norm(:,1),GFP_norm(:,1));
        [PscnI,tblI,statsI]=kruskalwallis([GFAP_norm(:,1),GFP_norm(:,1)]);
         [HscnI,PscnI]=kstest2(GFAP_norm(:,1),GFP_norm(:,1));
    end
end

 

%% choose data to plot:
if GFAP_OVX
    data_to_plot=per_positive_mean;
    %
    %
    % clear h
    % figure
    % subplot(2,1,1)
    % h=bar([1:2],nanmean(data_to_plot(Esr1_F,2:3)));h.FaceColor=[0.5 0.5 0.5];hold on % Females
    % plot([1*ones(size(data_to_plot(Esr1_F,2:3),1),1) 2*ones(size(data_to_plot(Esr1_F,2:3),1),1)],data_to_plot(Esr1_F,2:3),'k*')
    % h=bar([4 5],nanmean(data_to_plot(Esr1_FOVX,2:3),1));h.FaceColor=[0.5 0.5 0.5];hold on % Females
    % plot([4*ones(length(data_to_plot(Esr1_FOVX,2:3)),1) 5*ones(length(data_to_plot(Esr1_FOVX,2:3)),1)],data_to_plot(Esr1_FOVX,2:3),'k*')
    % xticks([1 2 4 5])
    % ylabel('mean integrated fluorescence')
    % xticklabels({'GFAP' 'Esr1' 'GFAP' 'Esr1'})
    % xlabel('F vs OVX')
    % ylim([0 32])
    %
    % subplot(2,1,2)
    % h=bar([1:2],nanmean(data_to_plot(PR_F,2:3),1));h.FaceColor=[0.5 0.5 0.5];hold on % Females
    % plot([1*ones(size(data_to_plot(PR_F,2:3),1),1) 2*ones(size(data_to_plot(PR_F,2:3),1),1)],data_to_plot(PR_F,2:3),'k*')
    % h=bar([4:5],nanmean(data_to_plot(PR_FOVX,2:3),1));h.FaceColor=[0.5 0.5 0.5];hold on % Females
    % plot([4*ones(length(data_to_plot(PR_FOVX,2:3)),1) 5*ones(length(data_to_plot(PR_FOVX,2:3)),1)],data_to_plot(PR_FOVX,2:3),'k*')
    % xticks([1 2 4 5])
    % ylabel('mean integrated fluorescence')
    % xticklabels({'GFAP' 'PR' 'GFAP' 'PR'})
    % xlabel('F vs OVX')
    % ylim([0 32])
    %
    %
    
    F=find(genders==1);% female
    FOVX=find(genders==2);% female OVX
    afi=union(F,FOVX);% all female index
    % normalied to max from females only
    % first index is mouse, second is channel
    norm_data_to_plot=[data_to_plot(:,1)./max(data_to_plot(afi,1)) data_to_plot(:,2)./max(data_to_plot(afi,2)) data_to_plot(:,3)./max(data_to_plot(afi,3))];
    
%     norm_data_to_plot_Esr1=[data_to_plot(:,1,:)./max(data_to_plot(afi,1,:)) data_to_plot(:,2,:)./max(data_to_plot(afi,2,:)) data_to_plot(:,3,:)./max(data_to_plot(union(Esr1_F,Esr1_FOVX),3,:))];
%     norm_data_to_plot_PR=[data_to_plot(:,1,:)./max(data_to_plot(afi,1,:)) data_to_plot(:,2,:)./max(data_to_plot(afi,2,:)) data_to_plot(:,3,:)./max(data_to_plot(union(PR_F,PR_FOVX),3,:))];
    
    % normalied to max from females only
  %  norm_data_to_plot=[data_to_plot(:,1,:)./max(data_to_plot(afi,1,:)) data_to_plot(:,2,:)./max(data_to_plot(afi,2,:)) data_to_plot(:,3,:)./max(data_to_plot(afi,3,:))];
%     
%     norm_data_to_plot_Esr1=[data_to_plot(:,1,:)./max(data_to_plot(afi,1,:)) data_to_plot(:,2,:)./max(data_to_plot(afi,2,:)) data_to_plot(:,3,:)./max(data_to_plot(union(Esr1_F,Esr1_FOVX),3,:))];
%     norm_data_to_plot_PR=[data_to_plot(:,1,:)./max(data_to_plot(afi,1,:)) data_to_plot(:,2,:)./max(data_to_plot(afi,2,:)) data_to_plot(:,3,:)./max(data_to_plot(union(PR_F,PR_FOVX),3,:))];
    
    %% figure that shows just GFAP: F vs OVX, Esr1 and PR, seperatly
    
    those_mice(F)
    those_mice(FOVX)
    
    %% plot just GFAP, Female vs OVX
    
    clear h
    figure
   
    
    h=bar(1,nanmean(norm_data_to_plot(F,2),1));h.FaceColor=[0.5 0.5 0.5];hold on % Females
    plot(1*ones(size(norm_data_to_plot(F,2),1),1),norm_data_to_plot(F,2),'ko')
    h=bar(2,nanmean(norm_data_to_plot(FOVX,2),1));h.FaceColor=[0.5 0.5 0.5];hold on % OVX
    plot(2*ones(length(norm_data_to_plot(FOVX,2)),1) ,norm_data_to_plot(FOVX,2),'ko')
    xticks([1 2])
    ylabel('mean integrated fluorescence')
    xticklabels({'F' 'OVX'})
    xlabel('GFAP')
    ylim([0 1.1])
   
    h2(1)=lillietest(norm_data_to_plot(F,2));
    h2(2)=lillietest(norm_data_to_plot(FOVX,2));
    if sum(h2)==2
        [HscnG,PscnG]=ttest2(norm_data_to_plot(F,2),norm_data_to_plot(FOVX,2));
    else
        %    [PscnG,HscnG]=ranksum(norm_data_to_plot(F,2),norm_data_to_plot(FOVX,2));
        [PscnG,tbl,stats]=kruskalwallis([norm_data_to_plot(F,2)',norm_data_to_plot(FOVX,2)'],[ones(1,length(F)),2*ones(1,length(FOVX))],'off');
    end
    
    
    % now Esr1
    
    
    %% GFAP, PR and ESr1 , female vs OVX
    
%     
%     clear h
%     figure
%     subplot(1,3,1)
%     h=bar(1,nanmean(norm_data_to_plot(F,2),1));h.FaceColor=[0.5 0.5 0.5];hold on % Females
%     plot(1*ones(size(norm_data_to_plot(F,2),1),1),norm_data_to_plot(F,2),'ko')
%     h=bar(2,nanmean(norm_data_to_plot(FOVX,2),1));h.FaceColor=[0.5 0.5 0.5];hold on % OVX
%     plot(2*ones(length(norm_data_to_plot(FOVX,2)),1) ,norm_data_to_plot(FOVX,2),'ko')
%     xticks([1 2])
%     ylabel('mean integrated fluorescence')
%     xticklabels({'F' 'OVX'})
%     xlabel('GFAP')
%     ylim([0 1.1])
%     [H,P]=ttest2(norm_data_to_plot(F,2),norm_data_to_plot(FOVX,2));
%     % now Esr1
%     clear h
%     ci=3;
%     subplot(1,3,2)
%     h=bar(1,nanmean(norm_data_to_plot_Esr1(Esr1_F,ci),1));h.FaceColor=[0.5 0.5 0.5];hold on % Females
%     plot(1*ones(size(norm_data_to_plot_Esr1(Esr1_F,ci),1),1),norm_data_to_plot_Esr1(Esr1_F,ci),'ko')
%     h=bar(2,nanmean(norm_data_to_plot_Esr1(Esr1_FOVX,ci),1));h.FaceColor=[0.5 0.5 0.5];hold on % Females
%     plot(2*ones(length(norm_data_to_plot_Esr1(Esr1_FOVX,ci)),1) ,norm_data_to_plot_Esr1(Esr1_FOVX,ci),'ko')
%     xticks([1 2])
%     ylabel('mean integrated fluorescence')
%     xticklabels({'F' 'OVX'})
%     xlabel('Esr1')
%     ylim([0 1.1])
%     [He,Pe]=ttest2(norm_data_to_plot_Esr1(Esr1_F,ci),norm_data_to_plot_Esr1(Esr1_FOVX,ci));
%     % now PR
%     clear h
%     
%     subplot(1,3,3)
%     h=bar(1,nanmean(norm_data_to_plot_PR(PR_F,ci),1));h.FaceColor=[0.5 0.5 0.5];hold on % Females
%     plot(1*ones(size(norm_data_to_plot_PR(PR_F,ci),1),1),norm_data_to_plot_PR(PR_F,ci),'ko')
%     h=bar(2,nanmean(norm_data_to_plot_PR(PR_FOVX,ci),1));h.FaceColor=[0.5 0.5 0.5];hold on % Females
%     plot(2*ones(length(norm_data_to_plot_PR(PR_FOVX,ci)),1) ,norm_data_to_plot_PR(PR_FOVX,ci),'ko')
%     xticks([1 2])
%     ylabel('mean integrated fluorescence')
%     xticklabels({'F' 'OVX'})
%     xlabel('PR')
%     ylim([0 1.1])
%     
%     
%     h(1)=lillietest(norm_data_to_plot_PR(PR_F,2));
%     h(2)=lillietest(norm_data_to_plot_PR(PR_FOVX,2));
%     if sum(h)==2
%         [HscnG,PscnG]=ttest2(norm_data_to_plot_PR(PR_F,2),norm_data_to_plot_PR(PR_FOVX,2));
%     else
%         %    [PscnG,HscnG]=ranksum(norm_data_to_plot_PR(PR_F,2),norm_data_to_plot_PR(PR_FOVX,2));
%         [PscnG,tbl,stats]=kruskalwallis([norm_data_to_plot_PR(PR_F,2),norm_data_to_plot_PR(PR_FOVX,2)],'off');
%     end
%     
%     
%     
%     [Hp,Pp]=ttest2(norm_data_to_plot_PR(PR_F,ci),norm_data_to_plot_PR(PR_FOVX,ci));
%     1
    




end








