function [results] = run_get_fluor_suite2P(exp_ID)
%% %% this function is used after get_data_suite2P, and inspection of cells that were also stained using LGS
% this function if for 'behavior'
% for Zstack use 'run_get_Zsatck_fluor_suite2P'
close all
clear norm_all_F all_F cell_ind F exp
NUM=[];TXT=[];RAW=[];
FIG=0;
SAVE_MAT=0; % for synchrony analysis in Python

show_cells=0;

if nargin==0
   % exp_ID='2170' ;  channel='647'% decreased rates
    % exp_ID='2174'; channel='647' % increase rates
   % exp_ID='Str22'; % no change/small decrease
    
    % exp_ID='Str25';% increase rates
    %  exp_ID='Str35';% increase rates
    %  exp_ID='Str38';% increase rates
    exp_ID='Str39'; channel='647'% increase rates 
    % exp_ID='Str39';channel='555'
    
    % exp_ID='Str24';% no change/ increased first 5 minutes 
   %  exp_ID='Str28';% no change
  % exp_ID='C57_1N';% increased 
    
end
%mouse_info.session_style='individual';
mouse_info.session_style='combined'; mouse_info.exp_per_sess=2; % suite2P was done on combined tiffs

this_path='C:\Users\anatk\Documents\Light_sectioning\';

switch exp_ID
    case '2170'
        mouse_info.mouse='WT2170';;mouse_info.M_num=2170;
        mouse_info.suite2P_mode='sparse';
        mouse_info.Z=960;
        exp='behavior'
        all_sess=[2:8];
        % all_sess=[1:2];
        
        for si=1:3; sess_names{si}=['saline ' num2str(si+1)];end
        for si=4:7; sess_names{si}=['cocaine ' num2str(si-3)];end
        %%%%%%%%/sess2 sess3 sess4 5     6     7     8
        %%cell
        switch mouse_info.session_style
            case 'individual'
                
                switch mouse_info.suite2P_mode
                    case 'non_sparse'
                        sess_stained_cells =[ 195   66    74    31    46    17    71;
                            197    159    25    39    31    22    46 ;
                            %30     23     23     49     37    nan    56;
                            142    135    nan    nan   nan     nan    86;
                            4       3      8      9     17      6    13;
                            nan     61     nan     nan   19    16    24;
                            70    121    170    127    35    nan    nan;
                            15   nan    14    14    11    nan   nan;
                            nan    nan    nan     98    52   nan   38]; % those cells are for non sparse mode
                        sess_nonstained_cells=[]; % not defined yet
                    case 'sparse'
                        sess_stained_cells =[ 195   66    74    31    46    17    71;
                            197    159    25    39    31    22    46 ;
                            % 30     23     23     49     37    nan    56;
                            142    135    nan    nan   nan     nan    86;
                            4       3      8      9    17      6    13;
                            nan     61    nan    nan   19    16    24;
                            70    121     170    127    35    nan    74;
                            15   nan      14     14    11    nan   nan;
                            nan    nan    nan     98    52   nan   38
                            119    134    144    140    78   65   nan];
                        sess_nonstained_cells=[17 11 6 5 5 20 30;
                            nan 9 16 6 3 3 2;
                            5 5 4 1 13 4 6;
                            nan 29 24 21 23 12 9;
                            10 8 12 13 8 5 8;
                            31 89 nan 43 57 35 19;
                            140 60 79 22 64 62 32 ];
                end
                
                % sess_stained_cells=1+[194,65,73,30,45,16,70;29,22,22,48,36,nan,55;3,2,7,8,16,5,12;45,16,57,102,89,10,25;
                %      69,120,nan,nan,34,nan,nan;1,1,1,37,32,12,10;28,36,36,38,30,21,45];%each raw is a cell
                
            case 'combined'
                sess_stained_cells =[ 1   2    3    5    6    17    21;
                    2    10    8    9    5    3    4 ];%%% this is not real
                sess_nonstained_cells=[2 5 6 21 31 4 5;
                    11 9 16 6 4 3 2];
                
        end
        %% load synchronosy data/Pegah Kassraian Fard 07/2020
        [NUM,TXT,~] = xlsread([this_path '\' mouse_info.mouse '\Synchrony_Analysis\Results_Synchrony2']);
    case '2174'
        mouse_info.mouse='WT2174'; exp='behavior' ;mouse_info.M_num=2174;
        mouse_info.Z=1100;
        all_sess=[2:8];
        all_sess=[1:2];
        mouse_info.suite2P_mode='sparse';
        for si=1:3; sess_names{si}=['saline ' num2str(si+1)];end
        for si=4:7; sess_names{si}=['cocaine ' num2str(si-3)];end
        %%%%%%%%/sess2 sess3 sess4 5     6     7     8
        %%cell
        switch mouse_info.session_style
            case 'individual'
                sess_stained_cells =[ 19 33 nan 52;
                    9 nan 30 87;
                    nan 17 30 46;
                    77 31 17 71;
                    45 25 12 nan;
                    87 31 16 13;
                    40 33 nan 61]';% transpose to 2170
                sess_nonstained_cells=[21 22 4 25 5 ;
                    nan 34 5 16 70;
                    12 29 25 21 9;
                    nan 50 22 25 4;
                    13 43 nan 32 1;
                    18 20 17 68 5;
                    16 63 nan 24 nan]';% transpose to 2170
                % sess=2; stained_cell_ind=[nan 19 33 66 52 ]; non_stained_cell_ind=[21 22 4 25 5 ];
                % sess=3; stained_cell_ind=[60 9 nan 30 87]; non_stained_cell_ind=[nan 34 5 16 70 ];
                %  sess=4;  stained_cell_ind=[72 nan 17 30 46]; non_stained_cell_ind=[12 29 25 21 9 ];
                %  sess=5;  stained_cell_ind=[35 77 31 17 71]; non_stained_cell_ind=[nan 50 22 25 4 ];
                %   sess=6; stained_cell_ind=[41 45 25 12 nan];non_stained_cell_ind=[13 43 nan 32 1 ];
                %   sess=7; stained_cell_ind=[41,87,31,16,13];non_stained_cell_ind=[18 20 17 68 5 ];
                % sess=8;  stained_cell_ind=[38 40 33 nan 61];non_stained_cell_ind=[16 63 nan 24 nan ];
            case 'combined'
                sess_stained_cells =[ 1   2    3    5    6    17    21;
                    2    10    8    9    5    3    4 ];%%% this is not real
                sess_nonstained_cells=[2 5 6 21 31 4 5;
                    11 9 16 6 4 3 2];%%% this is not real
                
        end
        
    case 'Str22'
        mouse_info.mouse='Str22'; exp='behavior' ;mouse_info.M_num=22;
        mouse_info.Z=935;
        all_sess=[1:2];
        mouse_info.suite2P_mode='sparse';%%% check here
        mouse_info.exp_per_sess=2;
        for si=1; sess_names{si}=['saline ' num2str(si)];end
        for si=2; sess_names{si}=['cocaine ' num2str(si-1)];end
        
        sess_stained_cells =[1 2 3 4 5 6 ;
            1 2 3 4 5 6]';%
        sess_nonstained_cells=[4 5 6 7 8 9;
            4 5 6 7 8 9 ]';%
    case 'Str25'
        mouse_info.mouse='Str25'; exp='behavior' ;mouse_info.M_num=25;
        mouse_info.Z=692;
        all_sess=[1:2];
        
        
        % session_style='individual';
        for si=1; sess_names{si}=['saline ' num2str(si)];end
        for si=2; sess_names{si}=['cocaine ' num2str(si-1)];end
        % for si=3; sess_names{si}=['saline and cocaine ' num2str(si-2)];end % 3 is a combined session- easier with cell identification
        sess_stained_cells =[1 2 3 4 5 2 ;
            1 2 3 4 5 2]';%
        sess_nonstained_cells=[4 5 6 1 4 2;
            4 5 6 1 4 2 ]';%
        
    case 'Str38'
        mouse_info.mouse='Str38'; exp='behavior' ;mouse_info.M_num=38;
        mouse_info.Z=1065;
        all_sess=[1:2];
        
        for si=1; sess_names{si}=['saline ' num2str(si)];end
        for si=2; sess_names{si}=['cocaine ' num2str(si-1)];end
        sess_stained_cells =[1 2 3 4 5 6 ;
            1 2 3 4 5 6]';%
        sess_nonstained_cells=[4 5 6 7 1 2;
            4 5 6 7 1 2 ]';%
    case 'Str35'
        mouse_info.mouse='Str35'; exp='behavior' ;mouse_info.M_num=35;
        mouse_info.Z=930;
        all_sess=[1:2];
        
        for si=1; sess_names{si}=['saline ' num2str(si)];end
        for si=2; sess_names{si}=['cocaine ' num2str(si-1)];end
        sess_stained_cells =[10 5 6 43 2 13 ;
            10 5 6 43 2 13]';%
        sess_nonstained_cells=[10 5 6 43 2 13;
            10 5 6 43 2 13]';%
        
        
    case 'Str24'
        mouse_info.mouse='Str24'; exp='behavior' ;mouse_info.M_num=24;
        mouse_info.Z=730;
        all_sess=[1:2];
        %mouse_info.session_style='individual';
        mouse_info.session_style='combined'; mouse_info.exp_per_sess=2; % suite2P was done on combined tiffs
        mouse_info.suite2P_mode='';
        for si=1; sess_names{si}=['saline ' num2str(si)];end
        for si=2; sess_names{si}=['cocaine ' num2str(si-1)];end
        sess_stained_cells =[10 5 6 43 2 13 ;
            10 5 6 43 2 13]';% numbers are not true
        sess_nonstained_cells=[10 5 6 43 2 13;
            10 5 6 43 2 13]';% numbers are not true
        
        
    case 'Str28'
        mouse_info.mouse='Str28'; exp='behavior' ;mouse_info.M_num=24;
        mouse_info.Z=680;
        all_sess=[1:2];
        %mouse_info.session_style='individual';
        mouse_info.session_style='combined'; mouse_info.exp_per_sess=2; % suite2P was done on combined tiffs
        mouse_info.suite2P_mode='';
        for si=1; sess_names{si}=['saline ' num2str(si)];end
        for si=2; sess_names{si}=['cocaine ' num2str(si-1)];end
        sess_stained_cells =[10 5 6 43 2 13 ;
            10 5 6 43 2 13]';% numbers are not true
        sess_nonstained_cells=[10 5 6 43 2 13;
            10 5 6 43 2 13]';% numbers are not true
        
    case 'Str39'
        mouse_info.mouse='Str39_LGS'; exp='behavior' ;mouse_info.M_num=39;
        mouse_info.Z=680;
        all_sess=[1:2];
        %mouse_info.session_style='individual';
        mouse_info.session_style='combined'; mouse_info.exp_per_sess=2; % suite2P was done on combined tiffs
        mouse_info.suite2P_mode='';
        for si=1; sess_names{si}=['saline ' num2str(si)];end
        for si=2; sess_names{si}=['cocaine ' num2str(si-1)];end
        % here sessions were analyzied together, as one movie, so cells ID
        % should be identifical for both sessions
      %  sess_stained_cells =[14 97 3 9 36 29 82 63 88 77 39 61;14 97 3 9 36 29 82 63 88 77 39 61]'; 
       % sess_nonstained_cells=[38 102 12 69 17 80 34 15 20 22 68 18; 38 102 12 69 17 80 34 15 20 22 68 18]';% 
               
       switch channel
           case '647'
               sess_stained_cells=[3 9 36 39 29 82 87;3 9 36 39 29 82 87]';
           case '555'
               sess_stained_cells=[82 87  ;82 87 ]';
       end
         sess_nonstained_cells=[38 68 12 69 17 80 18 15;38 68 12 69 17 80 18 15]';

                  
    case 'C57_1N'
        mouse_info.mouse='C571N'; exp='behavior' ;mouse_info.M_num=1;
        mouse_info.Z=600;
        all_sess=[1:2];
        %mouse_info.session_style='individual';
        mouse_info.session_style='combined'; mouse_info.exp_per_sess=2; % suite2P was done on combined tiffs
        mouse_info.suite2P_mode='';
        for si=1; sess_names{si}=['saline ' num2str(si)];end
        for si=2; sess_names{si}=['cocaine ' num2str(si-1)];end
        sess_stained_cells =[10 5 6 43 2 13 ;
            10 5 6 43 2 13]';% numbers are not true
        sess_nonstained_cells=[10 5 6 43 2 13;
            10 5 6 43 2 13]';% numbers are not true
end

switch exp_ID
    case {'2170' '2174'} 
        if strcmp(mouse_info.session_style,'individual')
        %sess_to_plot=[1,4:7];label_names={'sal#1' 'coc#1' 'coc#2' 'coc#3' 'coc#4'}; xtick_loc=[1 9 17 25 33];
        sess_to_plot=[1:7];label_names={'sal#1' 'sal#2' 'sal#3' 'coc#1' 'coc#2' 'coc#3' 'coc#4'}; xtick_loc=[1 9 17 25 33 41 49 57];
        %sess_to_plot=[1,7];label_names={'sal#1' 'coc#4'}; xtick_loc=[1 9];
    %case {'Str22', 'Str25', 'Str35', 'Str38'};
        else
           sess_to_plot=[1,2];label_names={'sal#1' 'coc#1'}; xtick_loc=[1 9];
        end
    otherwise
        sess_to_plot=[1,2];label_names={'sal#1' 'coc#1'}; xtick_loc=[1 9];
      %  sess_to_plot=[3];label_names={'sal and coc #1'}; xtick_loc=[1];
end

cells_ID='stained'
switch cells_ID
    case 'stained'
        cells_to_check=sess_stained_cells;
    case 'non-stained'
        cells_to_check=sess_nonstained_cells;
end

Fs=4;% sampling frame rate
per=2;

for si=1:length(all_sess); sess_markers{si}=['o'];end
% cell indexes, by eye. added one due to python indexing to matlab indexing
sess_fig=0;
switch mouse_info.session_style
    case 'individual'
        for si=1:length(all_sess)
            this function plot only the 'cell' , but returns fluorscence of all
            [norm_all_spks{si}, all_spks{si}, norm_all_dF{si},all_dF{si},all_F{si},cell_ind{si},cell_score{si}] = get_fluor_suite2P(mouse_info,all_sess(si),exp,sess_fig,Fs);
            arrays_length(si)=size(norm_all_dF{si},2);
        end
        for si=1:length(all_sess)
            on(si)=1;
            off(si)=min(arrays_length);
        end
    case 'combined'
         for si=1:length(all_sess)
             [norm_all_spks{si}, all_spks{si}, norm_all_dF{si},all_dF{si},all_F{si},cell_ind{si},cell_score{si}] = get_fluor_suite2P(mouse_info,all_sess(si),exp,sess_fig,Fs);
             arrays_length(si)=size(norm_all_dF{si},2)/2;
         end
         on=[1 min(arrays_length)] ;
         off=[1+min(arrays_length) 2*min(arrays_length)];

end
%% z-score  and set all same length
clear all_zscored_dF all_zscored_spks 
for si=1:length(all_sess)
    clear this_zscore_F norm_this_zscore_F this_zscore_spks 

    norm_all_spks{si}=norm_all_spks{si}(:,on(si):off(si));% set norm all spks to be the same size
    all_spks{si}=all_spks{si}(:,on(si):off(si));% set all spks to be the same size
    norm_all_dF{si}=norm_all_dF{si}(:,on(si):off(si));% set all F to be the same size
    all_dF{si}=all_dF{si}(:,on(si):off(si));% set all F to be the same size
   % all_F{si}=all_F{si}(:,on(si):off(si));% set all F to be the same size
    %this_dF=all_F{si};
    this_dF=all_dF{si};
    
    this_spks=all_spks{si};
    for ci=1:size(this_dF,1)
        clear X Y
        X=zscore(this_dF(ci,:));
        Y = prctile(X,per);
        this_zscore_dF(ci,:)=X-Y;
        norm_this_zscore_dF(ci,:)=this_zscore_dF(ci,:)./max(this_zscore_dF(ci,:));
        
        Xs=zscore(this_spks(ci,:));
        Ys=prctile(Xs,per);
        this_zscore_spks(ci,:)=Xs-Ys;
    end
    all_zscored_dF{si}=this_zscore_dF;
    all_norm_zscored_dF{si}=norm_this_zscore_dF;
    all_zscored_spks{si}=this_zscore_spks;
end


%Get the norm fluorescence of the reelevant cells
switch mouse_info.session_style
    case 'combined'
        sess_to_check=2;
    otherwise
       sess_to_check= size(cells_to_check,2);
end
dF_zscored=[];
for ci=1:size(cells_to_check,1)
    for si=1:sess_to_check
        if ~isnan(cells_to_check(ci,si))
            spks_zscored{ci}(si,:)=all_zscored_spks{si}(cells_to_check(ci,si),:);
            dF{ci}(si,:)=all_dF{si}(cells_to_check(ci,si),:);
%             F{ci}(si,:)=all_F{si}(cells_to_check(ci,si),:);
          
            dF_zscored{ci}(si,:)=all_zscored_dF{si}(cells_to_check(ci,si),:);
            % norm
            spks_norm{ci}(si,:)=norm_all_spks{si}(cells_to_check(ci,si),:);
            dF_norm{ci}(si,:)=norm_all_dF{si}(cells_to_check(ci,si),:);
%            F_norm{ci}(si,:)=norm_all_F{si}(cells_to_check(ci,si),:);
            dF_norm_zscored{ci}(si,:)=all_norm_zscored_dF{si}(cells_to_check(ci,si),:);
 %           F_norm_zscored{ci}(si,:)=all_norm_zscored_F{si}(cells_to_check(ci,si),:);
        else
            spks_zscored{ci}(si,:)=nan(1,min(arrays_length));
            dF{ci}(si,:)=nan(1,min(arrays_length));
         %   F{ci}(si,:)=nan(1,min(arrays_length));
            dF_zscored{ci}(si,:)=nan(1,min(arrays_length));
            % norm
            spks_norm{ci}(si,:)=nan(1,min(arrays_length));
            dF_norm{ci}(si,:)=nan(1,min(arrays_length));
            dF_norm_zscored{ci}(si,:)=nan(1,min(arrays_length));
        end
    end
end


  if SAVE_MAT
cd(['C:\Users\anatk\Documents\Light_sectioning\' mouse_info.mouse '\Synchrony_Analysis'])

      save ('sess_stained_cells.mat','sess_stained_cells')
      save('sess_nonstained_cells.mat','sess_nonstained_cells')
      save('all_norm_zscored_dF.mat','all_norm_zscored_dF')
      save('cell_ind.mat','cell_ind')
      mkdir('Data_files_code_v2')
      cd ('Data_files_code_v2')
      for si=1:size(sess_stained_cells,2)% by session
          k=0;m=0;
          for ci=1:size(sess_nonstained_cells,1)
              if ~isnan(sess_nonstained_cells(ci,si))
                  k=k+1;
                  this_nonst(k,:)=all_norm_zscored_dF{si}(sess_nonstained_cells(ci,si),:);
              end
              if ~isnan(sess_stained_cells(ci,si))
                  m=m+1;
                  this_st(k,:)=all_norm_zscored_dF{si}(sess_stained_cells(ci,si),:);
              end
          end
          this_cells=all_norm_zscored_dF{si}(cell_ind{si},:);
          eval(['sess' num2str(si+1) '_nonst' '=this_nonst;']);
          eval(['sess' num2str(si+1) '_st' '=this_st;']);
          eval(['sess' num2str(si+1) '=this_cells;']);
        
          
          save (['sess' num2str(si+1) '_nonst.mat'],['sess' num2str(si+1) '_nonst'])
          save (['sess' num2str(si+1) '_st.mat'],['sess' num2str(si+1) '_st'])
          save(['sess' num2str(si+1) '.mat'],['sess' num2str(si+1)])
      end
  end
 %% at this point Pegah is doing the syn analysis with python for 2170
 % this reads the results:
 if ~isempty(NUM)
     Syn_Matrix_all(:,1)=NUM([1:3 ,6:9],1);% all cells values, all 20 minutes
     Syn_Matrix_all(:,2)=NUM([1:3 ,6:9],2);% all cells SEM, all 20 minutes
     Syn_Matrix_all(:,3)=NUM([1:3 ,6:9],5);% all cells values,first  10 minutes
     Syn_Matrix_all(:,4)=NUM([1:3 ,6:9],6);% all cells SEM, first 10 minutes
     Syn_Matrix_all(:,5)=NUM([1:3 ,6:9],9);% all cells values,last  10 minutes
     Syn_Matrix_all(:,6)=NUM([1:3 ,6:9],10);% all cells SEM, last 10 minutes
     %% load ARC
     ARC_raws=[18,21,24,28,31,34,37];
     Syn_Matrix_ARC(:,1)=NUM(ARC_raws,1);% ARC cells values, all 20 minutes
     Syn_Matrix_ARC(:,2)=NUM(ARC_raws,2);% ARC cells SEM, all 20 minutes
     Syn_Matrix_ARC(:,3)=NUM(ARC_raws,5);% ARC cells values,first  10 minutes
     Syn_Matrix_ARC(:,4)=NUM(ARC_raws,6);% ARC cells SEM, first 10 minutes
     Syn_Matrix_ARC(:,5)=NUM(ARC_raws,9);% ARC cells values,last  10 minutes
     Syn_Matrix_ARC(:,6)=NUM(ARC_raws,10);% ARC cells SEM, last 10 minutes
     %% load non-ARC
     nonARC_raws=1+[18,21,24,28,31,34,37];
     Syn_Matrix_nonARC(:,1)=NUM(nonARC_raws,1);% nonARC cells values, all 20 minutes
     Syn_Matrix_nonARC(:,2)=NUM(nonARC_raws,2);% nonARC cells SEM, all 20 minutes
     Syn_Matrix_nonARC(:,3)=NUM(nonARC_raws,5);% nonARC cells values,first  10 minutes
     Syn_Matrix_nonARC(:,4)=NUM(nonARC_raws,6);% nonARC cells SEM, first 10 minutes
     Syn_Matrix_nonARC(:,5)=NUM(nonARC_raws,9);% nonARC cells values,last  10 minutes
     Syn_Matrix_nonARC(:,6)=NUM(nonARC_raws,10);% nonARC cells SEM, last 10 minutes
 end
%% plot synch analysis

if ~isempty(NUM)
    %% all cells
    fh=figure;
    % 20 minutes
    y=Syn_Matrix_all(sess_to_plot,1);sem=Syn_Matrix_all(sess_to_plot,2);
    x=[1:4:length(sess_to_plot)*4];
    bar(x,y,'FaceColor','w','EdgeColor','k','BarWidth',0.2,'LineWidth',1.5); hold on;
    er=errorbar(x,y,sem,sem); er.Color=[0 0 0]; er.LineStyle='none'; hold on;
    % first 10 minutes
    y=Syn_Matrix_all(sess_to_plot,3);sem=Syn_Matrix_all(sess_to_plot,4);
    x=[2:4:length(sess_to_plot)*4];
    bar(x,y,'FaceColor','w','EdgeColor',[0.2 0.2 0.2],'BarWidth',0.2,'LineWidth',1.5); hold on;
    er=errorbar(x,y,sem,sem); er.Color=[0 0 0]; er.LineStyle='none';  hold on;
    % last 10 minutes
    y=Syn_Matrix_all(sess_to_plot,5);sem=Syn_Matrix_all(sess_to_plot,6);
    x=[3:4:length(sess_to_plot)*4];
    bar(x,y,'FaceColor','w','EdgeColor',[0.4 0.4 0.4],'BarWidth',0.2,'LineWidth',1.5); hold on;
    er=errorbar(x,y,sem,sem); er.Color=[0 0 0]; er.LineStyle='none'; hold on;
    ylabel('Synchrony measure')
    title('all cells')
    set(gca,'XTick',[1 2 3 5 6 7 9 10 11 13 14 15 17 18 19])
    set(gca,'XTickLabel',{'20' 'F10' 'L10' '20' 'F10' 'L10' '20' 'F10' 'L10' '20' 'F10' 'L10' '20' 'F10' 'L10'})
    
    
    %% ARC vs nonARC 20minutes (full time)
    figure
    % 20 minutes ARC
    y=Syn_Matrix_ARC(sess_to_plot,1);sem=Syn_Matrix_ARC(sess_to_plot,2);
    x=[1:8:length(sess_to_plot)*8];
    bar(x,y,'FaceColor','w','EdgeColor',[0.5 0 0.5],'BarWidth',0.1,'LineWidth',1.5); hold on;
    er=errorbar(x,y,sem,sem); er.Color=[0 0 0]; er.LineStyle='none'; hold on;
    % 20 minutes nonARC
    y=Syn_Matrix_nonARC(sess_to_plot,1);sem=Syn_Matrix_nonARC(sess_to_plot,2);
    x=[2:8:length(sess_to_plot)*8];
    bar(x,y,'FaceColor','w','EdgeColor',[1 0.5 0],'BarWidth',0.1,'LineWidth',1.5); hold on;
    er=errorbar(x,y,sem,sem); er.Color=[0 0 0]; er.LineStyle='none'; hold on;
    ylabel('Synchrony measure')
    title('ARC vs nonARC')
    set(gca,'XTick',xtick_loc)
    set(gca,'XTickLabel',label_names)
    ylim([0 0.43])
    
    %% ARC vs non-ARC
    figure
    % first 10 minutes ARC
    y=Syn_Matrix_ARC(sess_to_plot,3);sem=Syn_Matrix_ARC(sess_to_plot,4);
    x=[3:8:length(sess_to_plot)*8];
    bar(x,y,'FaceColor','w','EdgeColor',[0.5 0 0.5],'BarWidth',0.1,'LineWidth',1.5); hold on;
    er=errorbar(x,y,sem,sem); er.Color=[0 0 0]; er.LineStyle='none';  hold on;
    % first 10 minutes nonARC
    y=Syn_Matrix_nonARC(sess_to_plot,3);sem=Syn_Matrix_nonARC(sess_to_plot,4);
    x=[4:8:length(sess_to_plot)*8];
    bar(x,y,'FaceColor','w','EdgeColor',[1 0.5 0],'BarWidth',0.1,'LineWidth',1.5); hold on;
    er=errorbar(x,y,sem,sem); er.Color=[0 0 0]; er.LineStyle='none'; hold on;
    % last 10 minutes
    y=Syn_Matrix_ARC(sess_to_plot,5);sem=Syn_Matrix_ARC(sess_to_plot,6);
    x=[5:8:length(sess_to_plot)*8];
    bar(x,y,'FaceColor','w','EdgeColor',[0.5 0 0.5],'BarWidth',0.1,'LineWidth',1.5); hold on;
    er=errorbar(x,y,sem,sem); er.Color=[0 0 0]; er.LineStyle='none'; hold on;
    % last 10 minutes nonARC
    y=Syn_Matrix_nonARC(sess_to_plot,5);sem=Syn_Matrix_nonARC(sess_to_plot,6);
    x=[6:8:length(sess_to_plot)*8];
    bar(x,y,'FaceColor','w','EdgeColor',[1 0.5 0],'BarWidth',0.1,'LineWidth',1.5); hold on;
    er=errorbar(x,y,sem,sem); er.Color=[0 0 0]; er.LineStyle='none'; hold on;
    
    ylabel('Synchrony measure')
    title('ARC vs nonARC F10minutes vs L10minutes')
    set(gca,'XTick',xtick_loc+3)
    set(gca,'XTickLabel',label_names)
    ylim([0 0.43])
end


%% plot dF
DATA='dF zscored'
%DATA='spikes'
%DATA='dF'
data_to_plot=[];
norm_to_plot=[];
switch DATA
    case 'spikes'
        data_to_calculate_events=all_zscored_spks;% all cells
        std_spikes=std2([all_zscored_spks{1} all_zscored_spks{2}]);
        data_to_plot=all_zscored_spks;
        norm_to_plot=all_zscored_spks;
    case 'dF'
        data_to_calculate_events=all_dF;% all cells
        data_to_plot=dF;
        norm_to_plot=dF_norm;
    case 'dF zscored' % zscored and perctile
        data_to_calculate_events=all_zscored_dF;% all cells
        std_spikes=std2([all_zscored_dF{1} all_zscored_dF{2}]);
        %data_to_calculate_events=all_zscored_spks; 
        if ~isempty(dF_zscored)
            data_to_plot=dF_zscored;
            norm_to_plot=dF_norm_zscored;
        end
end
 
t_array=[0:60/Fs:60*(min(arrays_length))/Fs]/3600;% time array in min
factor=1.5;


if ~isempty(cells_to_check) && show_cells
    %plot each cell by session
    for ci=1:size(cells_to_check,1)
        figure;
        subplot(2,2,1)
        for i=1:size(data_to_plot{ci},1)
            plot(t_array,data_to_plot{ci}(i,:)+i)
            hold on
        end
        xlabel('Time (min)')
        title([mouse_info.mouse ' ,cell ' num2str(cells_to_check(ci,1)) ' ' DATA ' ' cells_ID])
        
        subplot(2,2,3)
        imagesc(data_to_plot{ci})
        set(gca,'YDir','normal')
        
        subplot(2,2,2)
        for i=1:size(norm_to_plot{ci},1)
            plot(t_array,norm_to_plot{ci}(i,:)+i)
            hold on
        end
        xlabel('Time (min)')
        title([mouse_info.mouse ' ,cell '  num2str(cells_to_check(ci,1)) ' norm ' cells_ID])
        
        subplot(2,2,4)
        imagesc(norm_to_plot{ci})
        set(gca,'YDir','normal')
    end
    
    %COMPARE to cells that are found to be identical in all sessions using
    
    
    % check specific cell
    figure
    
    si=1;
    this_cell_ind=[1:size(cells_to_check,1)];
    for mi=1:length(this_cell_ind)
        subplot(length(this_cell_ind),1,mi)
        ci=cells_to_check(this_cell_ind(mi),si);
        if ~isnan(ci)
            peak_thresh=mean(data_to_calculate_events{si}(ci,:))+factor*std(data_to_calculate_events{si}(ci,:));
            findpeaks(data_to_calculate_events{si}(ci,:),t_array,'Annotate','extents','MinPeakProminence',double(peak_thresh));
            title([mouse_info.mouse ' ,sess ' num2str(si)+1 ' cell ' num2str(cells_to_check(mi,1)) ' ' cells_ID ' norm F'])
        end
        
    end
end
% find all cells peaks
for si=1:length(data_to_calculate_events)  
    for ci=1:size(data_to_calculate_events{si},1)
     %  figure
        %for ci=1:10
        %  subplot(10,1,ci)
     %    findpeaks(all_zscored_F{si}(ci,:),t_array,'Annotate','extents','MinPeakProminence',peak_thresh);
      %peak_thresh=mean(data_to_calculate_events{si}(ci,:))+factor*std(data_to_calculate_events{si}(ci,:));
      peak_thresh=mean(data_to_calculate_events{si}(ci,:))+factor*std_spikes;
      hL=round(0.5*length(data_to_calculate_events{si}(ci,:)));
      
        [allpks{si}{ci},alllocs{si}{ci},allw{si}{ci},allp{si}{ci}]=findpeaks(data_to_calculate_events{si}(ci,:),t_array,'Annotate','extents','MinPeakProminence',peak_thresh);
        [first10_allpks{si}{ci},first10_alllocs{si}{ci},first10_allw{si}{ci},first10_allp{si}{ci}]=findpeaks(data_to_calculate_events{si}(ci,1:hL),t_array(1:hL),'Annotate','extents','MinPeakProminence',peak_thresh);
        [last10_allpks{si}{ci},last10_alllocs{si}{ci},last10_allw{si}{ci},last10_allp{si}{ci}]=findpeaks(data_to_calculate_events{si}(ci,hL:end),t_array(hL:end),'Annotate','extents','MinPeakProminence',peak_thresh);

      %  is_burst(alllocs{si}{ci})
      %  hold on
      switch DATA
          case 'spikes'
              L=length(data_to_calculate_events{si}(ci,:));
%               F=factor*std(data_to_calculate_events{si}(ci,:));
%               F=0;
              F=factor*std_spikes;
              rates{si}(ci)=60*length(find(data_to_calculate_events{si}(ci,:)>F))/((min(arrays_length)-1)/Fs);% num peaks per minute
              first10_rates{si}(ci)=60*length(find(data_to_calculate_events{si}(ci,1:L/2)>F))/((min(arrays_length)-1)/Fs);% num peaks per minute
              last10_rates{si}(ci)=60*length(find(data_to_calculate_events{si}(ci,L/2+1:end)>F))/((min(arrays_length)-1)/Fs);% num peaks per minute
              
          case 'dF zscored'
              rates{si}(ci)=60*length(allpks{si}{ci})/((min(arrays_length)-1)/Fs);% num peaks per minute
              first10_rates{si}(ci)=60*length(first10_allpks{si}{ci})/(0.5*(min(arrays_length)-1)/Fs);% num peaks per minute
              last10_rates{si}(ci)=60*length(last10_allpks{si}{ci})/(0.5*(min(arrays_length)-1)/Fs);% num peaks per minute
      end
        peak_proms_mean{si}(ci)=mean(allp{si}{ci}.*allw{si}{ci}); % peak area 
        first10_peak_proms_mean{si}(ci)=mean(first10_allp{si}{ci}.*first10_allw{si}{ci}); % peak area 
        last10_peak_proms_mean{si}(ci)=mean(last10_allp{si}{ci}.*last10_allw{si}{ci}); % peak area 
        
        peak_width_mean{si}(ci)=mean(allw{si}{ci});
        first10_peak_width_mean{si}(ci)=mean(first10_allw{si}{ci});
        last10_peak_width_mean{si}(ci)=mean(last10_allw{si}{ci});
    end    
end


% plot event properties by session 
all_trait={'rates','rates first10', 'rates last10'};

%all_trait={'rates', 'area','width','rates first10', 'area first10','width first10','rates last10' 'area last10','width last10'};
%all_trait={'rates', 'area','width'};%,'rates first10', 'area first10','width first10','rates last10' 'area last10','width last10'};

plot_hist=0;
if plot_hist
    %% plot histograms of event rates
    % rates of selected cells positive features
    for ti=1:length(all_trait)
        trait=all_trait{ti}
        switch trait
            case 'rates'
                trait_to_plot=rates;LABEL='rate';hist_vect=[0:0.05:0.7];YMAX=7; XMAX=0.7
            case 'area'
                trait_to_plot=peak_proms_mean; LABEL='area'; hist_vect=[0:0.5:7]; YMAX=7; XMAX=7;
            case 'width'
                trait_to_plot=peak_width_mean; LABEL='width'; hist_vect=[0:10:200]; YMAX=6; XMAX=200;
        end
        figure
        for si=1:length(all_zscored_dF)
            subplot(2,ceil(length(all_zscored_dF)/2),si)
            % all defined cells/ non ARC cells
            indn1=sess_nonstained_cells(:,si); indn2=indn1(~isnan(indn1));
            histogram(trait_to_plot{si}(indn2),hist_vect,'FaceColor',[0 1 0])
            % histogram(rates{si}(cell_ind{si}),hist_vect,'FaceColor',[0 0 1])
            %   histogram(data_to_plot{si}(cell_ind{si}),'FaceColor',[0 0 1])
            hold on
            % just the ARC positive
            ind1=sess_stained_cells(:,si); ind2=ind1(~isnan(ind1));
            histogram(trait_to_plot{si}(ind2),hist_vect,'FaceColor',[0.5 0.5 0.5])
            % histogram(data_to_plot{si}(ind2),'FaceColor',[0 1 0])
            hold on
            % mean of all defined cells
            %line([mean(rates{si}(cell_ind{si})),mean(rates{si}(cell_ind{si}))],[0,YMAX],'Color',[0 0.1 1]);
            line([mean(trait_to_plot{si}(indn2)),mean(trait_to_plot{si}(indn2))],[0,YMAX],'Color',[0 1 0]);
            hold on
            % mean of ARC cells
            line([mean(trait_to_plot{si}(ind2)),mean(trait_to_plot{si}(ind2))],[0,YMAX],'Color',[0.5 0.5 0.5])
            title(['sess ' num2str(si+1)]); ylabel([ 'event ' LABEL ', all cells']); xlim([0 XMAX]);ylim([0 YMAX]);
            %  legend('non ARC');
        end
        
        % rates of ALL cell positive features
        figure
        for si=1:length(all_zscored_dF)
            subplot(2,ceil(length(all_zscored_dF)/2),si)
            % all defined cells/ non ARC cells
            histogram(trait_to_plot{si}(cell_ind{si}),hist_vect,'FaceColor',[0 1 0]); hold on
            %  histogram(rates{si}(cell_ind{si}),'FaceColor',[0 0 1])
            % mean of all defined cells
            line([mean(trait_to_plot{si}(cell_ind{si})),mean(trait_to_plot{si}(cell_ind{si}))],[0,35],'Color',[0 1 0.1]);hold on
            title(['sess ' num2str(si+1)]); ylabel([ 'event ' LABEL ', all cells']); xlim([0 XMAX]);ylim([0 YMAX]);
            xlim([0 XMAX]);ylim([0 YMAX]);
            %  legend('non ARC');
        end
    end
end
%%% statistics
switch exp_ID
    case {'2170','2174'}
        if ~isempty(sess_stained_cells)
            sess_ind1=1;
            %sess_ind2=7;
            sess_ind2=2;
            % index of cells, stained
            ind1a=sess_stained_cells(:,sess_ind1); ind2a=ind1a(~isnan(ind1a));% stained cells
            ind1b=sess_stained_cells(:,sess_ind2); ind2b=ind1b(~isnan(ind1b));% stained cells
            % index of cells, nonstained
            indn1a=sess_nonstained_cells(:,sess_ind1); indn2a=indn1a(~isnan(indn1a));% non stained cells
            indn1b=sess_nonstained_cells(:,sess_ind2); indn2b=indn1b(~isnan(indn1b));% non stained cells
            % define groups
            G1=first10_peak_width_mean{sess_ind1}(indn2a);
            G2=first10_peak_width_mean{sess_ind2}(indn2b);
            GN=[ones(1,length(G1)) 2*ones(1,length(G2))];
            % check KW test
            [pkw,~,~]=kruskalwallis([G1 G2],GN,'off');
            if pkw<0.05; hkw=1; else hkw=0; end
        end
end

%% plot trait using plotSpread
for ti=1:length(all_trait)
    trait=all_trait{ti};
    clear  sess_names_trait sess_markers_trait sess_markers_color
    switch trait
            case 'rates'
                trait_to_plot=rates; YMAX=60*0.25;%LABEL='rate';hist_vect=[0:0.05:0.7];YMAX=7; XMAX=0.7
            case 'area'
                trait_to_plot=peak_proms_mean; YMAX=1.7;%LABEL='prominance'; hist_vect=[0:0.5:7]; YMAX=7; XMAX=7;
            case 'width'
                trait_to_plot=peak_width_mean; YMAX=0.23;%LABEL='width'; hist_vect=[0:10:200]; YMAX=6; XMAX=200;
            case 'rates first10'
                trait_to_plot=first10_rates; YMAX=60*0.25;%LABEL='rate';hist_vect=[0:0.05:0.7];YMAX=7; XMAX=0.7
            case 'area first10'
                trait_to_plot=first10_peak_proms_mean; YMAX=1.7;%LABEL='prominance'; hist_vect=[0:0.5:7]; YMAX=7; XMAX=7;
            case 'width first10'
                trait_to_plot=first10_peak_width_mean; YMAX=0.23;%LABEL='width'; hist_vect=[0:10:200]; YMAX=6; XMAX=200;
           case 'rates last10'
                trait_to_plot=last10_rates; YMAX=60*0.25;%LABEL='width'; hist_vect=[0:10:200]; YMAX=6; XMAX=200;
          case 'area last10'
                trait_to_plot=last10_peak_proms_mean; YMAX=1.7;%LABEL='prominance'; hist_vect=[0:0.5:7]; YMAX=7; XMAX=7;
            case 'width last10'
                trait_to_plot=last10_peak_width_mean; YMAX=0.23;%LABEL='width'; hist_vect=[0:10:200]; YMAX=6; XMAX=200;


    end
    plot_info.trait=trait;
    plot_info.sess_names=sess_names;
    plot_info.sess_markers=sess_markers;
    plot_info.sess_to_plot=sess_to_plot;
    plot_info.YMAX=YMAX;
    plot_info.mouse_ID=mouse_info.mouse;
     plotSpread_cell_activity_trait_for_2P(trait_to_plot,plot_info,cell_ind,sess_stained_cells,sess_nonstained_cells)

end

%% check cells that had increased rates
switch mouse_info.mouse
    case 'Str39_LGS' % 2 session experiment, saline vs cocaine
rate_change=100*(rates{2}(cell_ind{1})-rates{1}(cell_ind{1}))./rates{1}(cell_ind{1});
increased=find(rate_change>=0);
decreased=find(rate_change<0);

increased_cell_ind=cell_ind{1}(increased); % cell_ind{1}==cell_ind{2} because the movies were analyzied in combined mode
decreased_cell_ind=cell_ind{1}(decreased);

stained_decreased=intersect(sess_stained_cells(:,1),decreased_cell_ind);
stained_increased=intersect(sess_stained_cells(:,1),increased_cell_ind);

[ii,inc_ind]=intersect(cell_ind{1},stained_increased);
[di,dec_ind]=intersect(cell_ind{1},stained_decreased);

find(ismember(increased,inc_ind)==0);

percent_rate_increase_all_nonARC=rate_change(increased(find(ismember(increased,inc_ind)==0)));
percent_rate_decrease_all_nonARC=rate_change(decreased(find(ismember(decreased,dec_ind)==0)));

percent_rate_increase_all_ARC=rate_change(inc_ind);
percent_rate_decrease_all_ARC=rate_change(dec_ind);

mean_percent_rate_increase_all_nonARC=mean(rate_change(increased(find(ismember(increased,inc_ind)==0))));
mean_percent_rate_decrease_all_nonARC=mean(rate_change(decreased(find(ismember(decreased,dec_ind)==0))));

mean_percent_rate_increase_all_ARC=mean(rate_change(inc_ind));
mean_percent_rate_decrease_all_ARC=mean(rate_change(dec_ind));

figure 
bar([1,2],[mean_percent_rate_increase_all_nonARC, mean_percent_rate_increase_all_ARC]); hold on;
bar([1,2],[mean_percent_rate_decrease_all_nonARC, mean_percent_rate_decrease_all_ARC])
ylabel('rate change (%)')


end


%% plot by mean cells activity by session 
if FIG
    for si=1:length(data_to_plot)
        data_per_session=[];
        for ci=1:size(cells_to_check,1)
            data_per_session=cat(1,data_per_session,data_to_plot{ci}(si,:));
        end
        mean_per_session(si,:)=mean(data_per_session,'omitnan');
    end
    
    figure
    for i=1:size(mean_per_session,1)
        plot(t_array,mean_per_session(i,:)+3*i-1)
        hold on
    end
    xlabel('Time (min)')
    ylabel('sessions')
    title([mouse_info.mouse ' ' cells_ID  ' cells, F norm'])
end
% plot cell 
if FIG
    figure
    for i=1:size(data_to_plot,2)
        plot(t_array,data_to_plot{i}+10*i-1); hold on
    end
end

plot_very_crowded_figure=1
all_trait2={'rates'};%,'area','width'};
if plot_very_crowded_figure
%% plot trait using bars, comparing first 10 minutes to last 10 minutes
for ti=1:length(all_trait2)
    trait=all_trait2{ti};
    clear data_to_plot sess_names_trait sess_markers_trait sess_markers_color trait_to_plot1
    switch trait
            case 'rates'
               trait_to_plot1=rates; YMAX=3.7;%LABEL='rate';hist_vect=[0:0.05:0.7];YMAX=7; XMAX=0.7
       
              %  trait_to_plot1=first10_rates; YMAX=60*0.05;%LABEL='rate';hist_vect=[0:0.05:0.7];YMAX=7; XMAX=0.7
                trait_to_plot2=rates; YMAX=3.7;%LABEL='rate';hist_vect=[0:0.05:0.7];YMAX=7; XMAX=0.7
            case 'area'
                  trait_to_plot1=peak_proms_mean; YMAX=0.8;%LABEL='prominance'; hist_vect=[0:0.5:7]; YMAX=7; XMAX=7;
           
             %   trait_to_plot1=first10_peak_proms_mean; YMAX=0.8;%LABEL='prominance'; hist_vect=[0:0.5:7]; YMAX=7; XMAX=7;
                trait_to_plot2=peak_proms_mean; YMAX=0.8;%LABEL='prominance'; hist_vect=[0:0.5:7]; YMAX=7; XMAX=7;
            case 'width'
                  trait_to_plot1=peak_width_mean; YMAX=0.23;%LABEL='width'; hist_vect=[0:10:200]; YMAX=6; XMAX=200;
         
                %trait_to_plot1=first10_peak_width_mean; YMAX=0.23;%LABEL='width'; hist_vect=[0:10:200]; YMAX=6; XMAX=200;
                trait_to_plot2=peak_width_mean; YMAX=0.23;%LABEL='width'; hist_vect=[0:10:200]; YMAX=6; XMAX=200;
    end
    iscell_trait_array1=[];
     iscell_trait_array2=[];
    isecell_trait_group=[];
 %   sess_to_plot=[1 3 4 7];% saline 2,4 coc 2,4
    % sess_to_plot=[1,4:7];% 
     sess_to_plot=[1 2];% 
    for si=1:length(sess_to_plot)
        sess_ind=sess_to_plot(si);
        iscell_trait1{si}=trait_to_plot1{sess_ind}(cell_ind{sess_ind});
        iscell_trait2{si}=trait_to_plot2{sess_ind}(cell_ind{sess_ind});
        iscell_trait_array1=[iscell_trait_array1 trait_to_plot1{sess_ind}(cell_ind{sess_ind})];
        iscell_trait_array2=[iscell_trait_array2 trait_to_plot2{sess_ind}(cell_ind{sess_ind})];
        isecell_trait_group=[isecell_trait_group sess_ind*ones(1,length(trait_to_plot1{sess_ind}(cell_ind{sess_ind})))];
        ind1=sess_stained_cells(:,sess_ind); ind2=ind1(~isnan(ind1));% stained cells
        % isStained_trait{si}=data_to_plot{si}(ind2);
        indn1=sess_nonstained_cells(:,sess_ind); indn2=indn1(~isnan(indn1));% non stained cells
        % isNotStained_trait{si}=data_to_plot{si}(indn2);
        specific_cells_trait1{si}=trait_to_plot1{sess_ind}(ind2);% ARC
        specific_cells_trait2{si}=trait_to_plot2{sess_ind}(indn2);% nonARC
        specific_cells_trait2{si}=trait_to_plot2{sess_ind}(indn2);% nonARC
%         non_specific_cells_trait1{si*2}=trait_to_plot1{sess_ind}(indn2);
%         non_specific_cells_trait2{si*2}=trait_to_plot2{sess_ind}(indn2);
        sess_names_trait{si}=[sess_names{sess_ind} ' ARC'];
        sess_names_trait2{si}=[sess_names{sess_ind} ' nonARC'];
        sess_markers_trait{si}='o';
        sess_markers_trait2{si}='*';
        sess_markers_color{si}=[0.5 0 0.5];
        sess_markers_color2{si}=[1 0.5 0];
    end
    
    for si=1:length(iscell_trait1)
        for si2=si:length(iscell_trait1)
            % first 10 minutes
            [h1(si,si2),p1(si,si2)]=kstest2(iscell_trait1{si},iscell_trait1{si2});
            [pkw1(si,si2),~,~]=kruskalwallis([iscell_trait1{si} iscell_trait1{si2}],[ones(1,length(iscell_trait1{si})) 2*ones(1,length(iscell_trait1{si2}))],'off');
            if pkw1(si,si2)<0.05; hkw1(si,si2)=1; else hkw1(si,si2)=0; end
             % last 10 minutes
             [h2(si,si2),p2(si,si2)]=kstest2(iscell_trait2{si},iscell_trait2{si2});
            [pkw2(si,si2),~,~]=kruskalwallis([iscell_trait2{si} iscell_trait2{si2}],[ones(1,length(iscell_trait2{si})) 2*ones(1,length(iscell_trait2{si2}))],'off');
            if pkw2(si,si2)<0.05; hkw2(si,si2)=1; else hkw2(si,si2)=0; end
  
        end
    end
    
    
    for si=1:length(specific_cells_trait1)
        for si2=si:length(specific_cells_trait1)
            % ARC
            [hsc1(si,si2),psc1(si,si2)]=kstest2(specific_cells_trait1{si},specific_cells_trait1{si2});
            [pkwsc1(si,si2),~,~]=kruskalwallis([specific_cells_trait1{si} specific_cells_trait1{si2}],[ones(1,length(specific_cells_trait1{si})) 2*ones(1,length(specific_cells_trait1{si2}))],'off');
            if pkwsc1(si,si2)<0.05; hkwsc1(si,si2)=1; else hkwsc1(si,si2)=0; end
            
                  % non ARC
            [hsc2(si,si2),psc2(si,si2)]=kstest2(specific_cells_trait2{si},specific_cells_trait2{si2});
            [pkwsc2(si,si2),~,~]=kruskalwallis([specific_cells_trait2{si} specific_cells_trait2{si2}],[ones(1,length(specific_cells_trait2{si})) 2*ones(1,length(specific_cells_trait2{si2}))],'off');
            if pkwsc2(si,si2)<0.05; hkwsc2(si,si2)=1; else hkwsc2(si,si2)=0; end

        end
    end
   % p = kruskalwallis(iscell_trait_array,isecell_trait_group);
   % mean of all cells
   for fi=1:length(iscell_trait1)
       mean_iscell_trait1(fi)=nanmean(iscell_trait1{fi});
       sem_iscell_trait1(fi)=std(iscell_trait1{fi})/sqrt(length(iscell_trait1{fi}));
        mean_iscell_trait2(fi)=nanmean(iscell_trait2{fi});
       sem_iscell_trait2(fi)=std(iscell_trait2{fi})/sqrt(length(iscell_trait2{fi}));
   end
   %% mean for ARC/non-ARC
   for fi=1:length(specific_cells_trait1)
       mean_specific_cell_trait1(fi)=nanmean(specific_cells_trait1{fi});%ARC
       sem_specific_cell_trait1(fi)=std(specific_cells_trait1{fi})/sqrt(length(specific_cells_trait1{fi}));
       mean_specific_cell_trait2(fi)=nanmean(specific_cells_trait2{fi});% nonARC
       sem_specific_cell_trait2(fi)=std(specific_cells_trait2{fi})/sqrt(length(specific_cells_trait2{fi}));
   end
    % show distributions of event RATES using plotSpread function or bar
    % plot
    figure
    %subplot(1,2,1)
    % plot the first 10 minutes
    x1=[1:length(mean_iscell_trait1)];
    bar(x1,mean_iscell_trait1,'FaceColor','w','EdgeColor',[0 1 0],'BarWidth',0.1,'LineWidth',1.5); hold on;
    er=errorbar(x1,mean_iscell_trait1,sem_iscell_trait1,sem_iscell_trait1); er.Color=[0 0 0]; er.LineStyle='none'; hold on;
    % plot the last 10 minutes
%     x2=[2:3:3*length(mean_iscell_trait2)];
%     bar(x2,mean_iscell_trait2,'FaceColor','w','EdgeColor',[0.1 0.8 0.2],'BarWidth',0.1,'LineWidth',1.5); hold on;
%     er=errorbar(x2,mean_iscell_trait2,sem_iscell_trait2,sem_iscell_trait2); er.Color=[0 0 0]; er.LineStyle='none'; hold on;

    % plot a line when is significant 
    step=0;
    for si=1:length(iscell_trait1)
        for si2=si:length(iscell_trait1)
            if hkw1(si,si2)
                step=step+0.01*YMAX;
               line([x1(si) x1(si2)],[YMAX-step YMAX-step]); hold on
            end
            
%             if hkw2(si,si2)
%                 step=step+0.01*YMAX;
%                line([x2(si) x2(si2)],[0.8*YMAX-step 0.8*YMAX-step]); hold on
%             end
        end
    end
   
  
   %  ylim([0 YMAX])
    %%% plot stained cells vs non stained
   %subplot(1,2,2)
%     plotSpread(specific_cells_trait, ...
%         'xNames', sess_names_trait, ...
%         'distributionMarkers', sess_markers_trait,'distributionColors',sess_markers_color);
%     hold on
   x2=[x1(2)+2:x1(2)+1+length(mean_specific_cell_trait1)];%ARC
    bar(x2,mean_specific_cell_trait1,'FaceColor','w','EdgeColor',[0 1 0],'BarWidth',0.1,'LineWidth',1.5); hold on;
    er=errorbar(x2,mean_specific_cell_trait1,sem_specific_cell_trait1,sem_specific_cell_trait1); er.Color=[0 0 0]; er.LineStyle='none'; hold on;
    % nonARC
    if strcmp(channel,'647')
        x3=[x2(2)+2:x2(2)+1+length(mean_specific_cell_trait2)];
        bar(x3,mean_specific_cell_trait2,'FaceColor','w','EdgeColor',[0.1 0.8 0.2],'BarWidth',0.1,'LineWidth',1.5); hold on;
        er=errorbar(x3,mean_specific_cell_trait2,sem_specific_cell_trait2,sem_specific_cell_trait2); er.Color=[0 0 0]; er.LineStyle='none'; hold on;
    end
      % plot a line when is significant 
    step=0;
    for si=1:length(mean_specific_cell_trait1)
        for si2=si:length(mean_specific_cell_trait1)
            if hkwsc1(si,si2)
                step=step+0.01*YMAX;
               line([x2(si) x2(si2)],[YMAX-step YMAX-step]); hold on
            end
            
%             if hkwsc2(si,si2)
%                 step=step+0.01*YMAX;
%                line([x2(si) x2(si2)],[0.8*YMAX-step 0.8*YMAX-step]); hold on
%             end
        end
    end
      ylabel(['peak event ' trait ])
     ylim([0 YMAX])
       title(['peak event ' trait ' all ' mouse_info.mouse  ' cahnnel ' channel])
 
end
end


%plot specific sess by order, choosing the session that had a significant
%different between stained to non stained
%optional output to plot: all_zscored_dF, norm_all_dF, all_norm_zscored_dF
sess=sess_to_plot;
sess=[1 2];
clear B I inds
%define the order of cells
 BY='rates'
 BY_sess=1;
inds=cell_ind{sess(2)};
[B,I]=sort(rates{sess(2)}(inds));
% inds without the ARC positive cells
Is=ismember(inds(I), sess_stained_cells(:,BY_sess));
A=inds(I); A=A(Is==0) ;
% combine with the ARC positive cells
myind=[sess_stained_cells(:,BY_sess)' A'];


figure
for si=1:length(sess)  
    subplot(1,max(sess)+1,si)
    % imagesc(norm_all_F{sess(si)}(cell_ind{sess(si)},:))
    clims = [0 1];
    imagesc(t_array,[0 size(I)] ,norm_all_dF{sess(si)}(flip(I),:),clims)
    %imagesc(t_array,[0 size(I)] ,all_norm_zscored_dF{sess(si)}(I,:),clims)
    colormap('hot')
    colorbar
    title([mouse_info.mouse ' ,sess ' num2str(sess(si)+1) ' dF map'])
    xlabel('Time (min)')
    ylabel ([ 'all cells (sorted by ' BY ')'])
end
switch mouse_info.session_style % add 
    case 'combined'
        
        full_t_array=[0:60/Fs:60*2*(min(arrays_length))/Fs]/3600;% time array in min
        norm_cell_dF=[norm_all_dF{1}(myind,:) norm_all_dF{2}(myind,:)];
        
        subplot(1,max(sess)+1,max(sess)+1)
       % figure
        for i=1:size(norm_cell_dF,1)
            plot(full_t_array,norm_cell_dF(i,:)+size(norm_cell_dF,1)-i)
            hold on
        end
        title(['norm dF'])
        ylim([0 size(norm_cell_dF,1)+1])
end

1


if FIG
clear B I inds
all_which_cells={'non_stained','stained'};

for i=1:length(all_which_cells)
    figure
    sort_by_sess=2;
    switch all_which_cells{i}
        case 'non_stained'
            inds=sess_nonstained_cells(:,sort_by_sess);
            %inds=inds(~isnan(inds));
        case 'stained'
            inds=sess_stained_cells(:,sort_by_sess);
            %inds=inds(~isnan(inds));
    end
    BY='rate';
    [B,I]=sort(rates{sort_by_sess}(inds(~isnan(inds))));
    for si=1:length(sess)
        
        which_cells=all_which_cells{i};
        switch which_cells
            case 'non_stained'
                inds=sess_nonstained_cells(:,sess(si)); 
            case 'stained'
                inds=sess_stained_cells(:,sess(si));
        end
        
        subplot(1,length(sess),si)
        % imagesc(norm_all_F{sess(si)}(cell_ind{sess(si)},:))
        tmpinds=inds(I); inds=tmpinds(~isnan(tmpinds));
        imagesc(t_array,[0 size(I)] ,all_norm_zscored_dF{sess(si)}(inds,:))
         switch which_cells; case 'non_stained'; colormap('autumn'); case 'stained'; colormap('pink'); end 
         switch which_cells; case 'non_stained'; colormap('hot'); case 'stained'; colormap('hot'); end 
        title([mouse_info.mouse ' ,sess ' num2str(sess(si)+1) ' norm F map'])
        xlabel('Time (min)')
        ylabel ([which_cells ' cells (sorted by ' BY ')'])
    end
end
end
% 
% for si=1:length(all_sess)
%     sess_ind=all_sess(si);
% 
%     results.rates{si}=rates{si}(cell_ind{si});
% 
%     results.peak_proms_mean{si}=peak_proms_mean{si}(cell_ind{si});
%     results.peak_width_mean{si}=peak_width_mean{si}(cell_ind{si});
%     results.first10_rates{si}=first10_rates{si}(cell_ind{si});
%     results.first10_peak_proms_mean{si}=first10_peak_proms_mean{si}(cell_ind{si});
%     results.first10_peak_width_mean{si}=first10_peak_width_mean{si}(cell_ind{si});
%     results.last10_rates{si}=last10_rates{si}(cell_ind{si});
%     
%     results.sess_stained_cells(si)=sess_stained_cells(:,si);
%     results.sess_nonstained_cells(si)=sess_nonstained_cell(:,si);
%     
% end
% 

clear results
results.cell_ind=cell_ind;
results.sess_stained_cells=sess_stained_cells;
results.sess_nonstained_cells=sess_nonstained_cells;
results.rates=rates;
results.peak_proms_mean=peak_proms_mean;
results.peak_width_mean=peak_width_mean;
results.first10_rates=first10_rates;
results.first10_peak_proms_mean=first10_peak_proms_mean;
results.first10_peak_width_mean=first10_peak_width_mean;
results.last10_rates=last10_rates;
results.last10_peak_proms_mean=last10_peak_proms_mean;
results.last10_peak_width_mean=last10_peak_width_mean;
results.all_trait=all_trait;

1

