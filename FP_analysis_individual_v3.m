function [all_dF] = FP_analysis_individual_v3(ID,side,date,Sess,Sname,Gender,Estrus,rig,MinPeakP,interval,bins,thresh)
% Anat Kahan, Cell Reports, 2021
experiment_type='VIPFP'
%experiment_type='VIPFP_Onset'
%experiment_type='iCreV'
clear files_name1 files1 all_dF
% LFcut = 4; % cut-off frequency of lowpass filter
% order = 4; % N-th order for butterworth filter
POWER=0;
%% 
%mypath='D:\Data_Glab\';
mypath='D:\Data_Glab_home_work\'; % external drive


FIGbinned=1;
FIG=0;
if nargin == 0
    
    FIG=1;
      
    ID='VIPGC103R'; side='R'; Gender='M'; rig='TDT';
    % ID='VIPGC106LL'; side='R'; Gender='M'; rig='SynTDT'; 
    %ID='VIPGC122R'; side='R'; Gender='F'; rig='SynTDT';
    % ID='VIPGC123L'; side='R'; Gender='F'; rig='SynTDT';
   % ID='VIPGC198L'; side='R'; Gender='F'; rig='SynTDT';
    
    
   % date='090318';
     date='042819';
  % date='051520'
   
    Sname='test6R1';%  
    Estrus='na';
       
end

if strfind(Sname,'test')
    par_file_name='ParametersFPdFOnly';
else
    par_file_name='ParametersFP';
end
[NUMpar,TXTpar,RAWpar]=xlsread([mypath 'fiberphotometry\' par_file_name '.xlsx']);
rel=NUMpar(1,7); % if peak analysis 'MinPeak' is relative to specific trial
interval=NUMpar(1,1);
bins=NUMpar(:,2)';%[2,5,10,20,30]; %time in minutes
thresh=NUMpar(1,4);%3;

% set more parameters
if strfind(Sname,'test') 
    Zscore=0; testYLIM=35;
    %NUMpar(1,6)=0;
elseif strfind(Sname,'Lumencore')
    %NUMpar(1,6)=0;
    Zscore=0; 
else
    Zscore=NUMpar(1,6); 
end
USE_INT_FOR_dF=0;
Smth=1;
Perc=NUMpar(1,5);
Lpass=1;
Hpass=1;
%MinPeakP=4;% peak analysis parameter
 MinPeakP=NUMpar(1,3); % peak analysis parameter 
MinPeakW=4;% in seconds. later that should be read from xls file

switch rig
    case 'TDT'
        path=[mypath 'fiberphotometry\TDT_FP\'];
    case 'SynTDT'
        path=[mypath 'fiberphotometry\SynTDT_FP\'];
end
    
files_name1=[ID ' ' Gender ' ' Estrus ' ' side 'fiber ' date  ' '  Sname]
files1=[ID '_' side 'fiber_' date '_' Sname];

% for SCN VIP experiments, in which TTL show when light goes off


%% read the df data
fullpath=[path files1];
load(fullpath);
all_dF= y;
if strcmp(files1,'VIPGC103R_Rfiber_042819_test1'); all_dF.data=[all_dF.data(2,:); all_dF.data(1,:)];end
is_onset.is=strmatch('SessOnset',Sname);
if is_onset.is
    files2=[ID '_' side 'fiber_' date '_Sess' Sname(10:end) ];
    is_onset.L=length(all_dF.t);
    if exist([path files2])
    load([path files2]);
    all_dF2=y;
    % combine Onset and dF to get a better treatmentment for baseline etc. 
    all_dF.data=[all_dF.data all_dF2.data];
    all_dF.t=[all_dF.t all_dF2.t];
    all_dF.dF=[all_dF.dF all_dF2.dF];
    all_dF.fit400=[all_dF.fit400 all_dF2.fit400];
    end
end

switch rig
    case 'TDT'
        fs=382; % TDT FP rig
    case 'SynTDT'
        fs = y.fs; %Syn TDT FP rig 
end

t1=fs*10; %skipping the first 10 seconds, the freq and intensity set
t1=fs*2; %skipping the first 2 seconds, the freq and intensity set
%all_dF.Data= all_dF.Data(:,t1:end);
switch files1
    case 'VIPGC128R_Rfiber_062419_test6R1'
         t1=fs*1;
          t2=length(all_dF.data);
    otherwise
        t2=length(all_dF.data);
end

%% finds light status
[light_array]=finds_light_status(files1,all_dF.TTL,Sname,t1,fs);

% define relevant time 08/06/2019
switch light_array.exp
    case {'test6R','test'}
        TRANGE=[0 300]; time_epoc=0;
        %isonset=0;
end

[dF,t,all_dF]=get_df_from_raw_data2(all_dF,fs,NUMpar,t1,t2,TRANGE,time_epoc,is_onset);
    
%% check
%figure;plot(t,dF);hold on; plot(all_dF.t,all_dF.dF);title(files_name1)
t_original=all_dF.t;
all_dF.t=t;
all_dF.dF=dF;


%% reduce sampling rate for data and fit
%reducesampling_f=0.05;
red_data(1,:)=interp1(t_original,all_dF.data(1,:),[1:0.05:t_original(end)]);
red_data(2,:)=interp1(t_original,all_dF.data(2,:),[1:0.05:t_original(end)]);
all_dF.data=red_data;
new_fs=length(all_dF.dF)/all_dF.t(end);
%all_dF.fit400=interp1(all_dF.t,all_dF.fit400,[1:0.05:all_dF.t(end)]);
all_dF.fit400_original=interp1(t_original,all_dF.fit400_original,[1:0.05:t_original(end)]);
peak_thresh=MinPeakP;
disp(['This sample peak thresh is: ' num2str(peak_thresh)]) 

%% find peaks using 'findpeaks' matlab function
% index 1 - all, 2- dark phase, 3 - light phase
state_str={'all' ;'dark phase'; 'light phase'};

%[pks,locs,w,p]=findpeaks(all_dF(k,i).dF,all_dF(k,i).t,'Annotate','extents','WidthReference','halfheight','MinPeakProminence',3);
this_Sname=Sname(1:end-1);
if str2num(this_Sname(end))>0
  this_Sname=this_Sname(1:end-2);  
end  
switch this_Sname
    case {'test6R','test'}
        [pks{1},locs{1},w{1},p{1}]=findpeaks(all_dF.dF,all_dF.t,'Annotate','extents','WidthReference','halfheight','MinPeakProminence',peak_thresh,'MinPeakWidth',MinPeakW);
        state_str={'all' };
        for li=1:length(light_array.light_on)
            this_light_ind=intersect(find(all_dF.t>light_array.light_on(li)),find(all_dF.t<light_array.light_off(li)));
            this_dark_ind=intersect(find(all_dF.t<light_array.light_on(li)),find(all_dF.t>light_array.light_on(li)-15)); % 5 sec before light on is baseline
        
            dF_ON(li)=sum(all_dF.dF(this_light_ind));
            dF_OFF(li)=sum(all_dF.dF(this_dark_ind));
            peak_amplitude(li)=max(all_dF.dF(this_light_ind))-max(all_dF.dF(this_dark_ind));
          %  figure; plot(all_dF.t(this_light_ind),all_dF.dF(this_light_ind));
        end       
end

for pi=1:length(state_str)
    num_picks(pi)=length(pks{pi});
    peaks_area(pi)=mean(w{pi}.*p{pi});
    peaks_height(pi)=mean(p{pi});
end
switch this_Sname
     case {'test6R','test'}
           peaks_area(pi)=nanmean(dF_ON)-nanmean(dF_OFF);
           peaks_height(pi)=nanmean(peak_amplitude);  
end

%%
session_length_sec(1)=max(all_dF.t);
switch this_Sname
    case {'Ses','Sess','DL','testR6','LDold','DLold'}
        session_length_sec(2)=TRANGE(2);
        session_length_sec(3)=TRANGE(2);
end

peak_analysis.state_str=state_str;
peak_analysis.num_picks=num_picks;
peak_analysis.session_length_sec=session_length_sec;
peak_analysis.peaks_area=peaks_area;
peak_analysis.peaks_height=peaks_height;
all_dF.light_array=light_array;
all_dF.peak_analysis=peak_analysis;
all_dF.bins=bins;
MinPeakPstr=num2str(MinPeakP); pind=strfind(MinPeakPstr,'.');
if abs(MinPeakP-round(MinPeakP))>0;MinPeakPstr=[MinPeakPstr(1:pind-1) 'p' MinPeakPstr(pind+1:end)];end
% now save
all_dF.Data=[];
all_dF.fit400=[];
all_dF.TTL=[];

if rel==1
    save([fullpath '_int_processed' MinPeakPstr 'rel_' num2str(interval) '_Perc' num2str(Perc) '_ZS' num2str(Zscore)],'all_dF')
elseif rel==0
     save([fullpath '_int_processed' MinPeakPstr '_' num2str(interval) '_Perc' num2str(Perc) '_ZS' num2str(Zscore)],'all_dF')
end  

if FIG
    %% dF/F figure
    switch this_Sname
        case 'test6R'
            XMAX=320;
   end
    
    figure
    subplot (2,1,1)
    findpeaks(all_dF.dF,all_dF.t,'Annotate','extents','MinPeakProminence',peak_thresh);    hold on
    ylim([-5,15]);
    xlim([0,XMAX]);
  
    title([files_name1 ' ,Min Peak=' num2str(MinPeakP)])
    xlabel('Time (sec)','FontSize',12);
    ylabel('dF (Z-score)','FontSize',12);
    
    subplot (2,1,2)
    plot(all_dF.t,all_dF.data(1,:));     hold on
    plot(all_dF.t,all_dF.data(2,:));     hold on
    xlim([0,XMAX]);
   
    xlabel('Time (sec)','FontSize',12);
    ylabel('raw data','FontSize',12);
    
    
end



1
end


%%%%%%%%%%





%% fit function, to remove ref from signal
function [dF] = fit_ref(data)
    B = data(1,:)';
    A = [data(2,:)' ones(length(data),1)];
    theta = A\B;
    fit400 = theta(1)*data(2,:)+theta(2);
    dF = 100*detrend((data(1,:)-fit400)./fit400);
end

function y = lowpass(x,cutoff,Fs,n)
% Butterworth Lowpass Filter (zero-phase distortion filter)
% This creates n-th order Butterworth lowpass filter and takes input
% signal and creates output, the filtered signal. 
%
% <Usage>
%
% y = lowpass(x,cutoff,Fs,n)
% 
% where x: the unfiltered, raw signal
%       cutoff: cut-off frequency
%       Fs: sampling rate
%       n: the order of Butterworth filter
%       y: the filtered signal
%
% <Example>
%
% y = lowpass(x,100,2000,4);
%
% Coded by Ryan Cho, Oct 21 2013

[b,a] = butter(n,cutoff/(Fs/2),'low');
y = filtfilt(b,a,double(x));
end

function y = highpass(x,order,cutoff,Fs)
    [z,p,k] = butter(order,cutoff/Fs,'high');
    y = filtfilt(z,p,double(x));
end
