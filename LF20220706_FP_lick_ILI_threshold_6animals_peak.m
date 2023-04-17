%% Script written by Lauren Faget - Hnasko lab - July 2021
% For Fiber photometry data analysis

% Make sure to have MATLAB scripts, MATLAB tools and TDT tanks in the same folder. 
% In 'Home' tab, click on 'Set Path', click on 'Add with subfolders' and
% select your folder.

% This script average data for a group of 8 animals and 4 epocs (TTL input types)
% If less than 8 animals and 4 epocs, comment unnecessary sections 
% comment : ctrl R
% uncomment : ctrl T

%% Useful commands:
%clear --> clear the workspace data
%clc --> clear the command window
%close all --> close all graphs open
%% 
clc
clear
%%
% Script can be run all at once by ckicking on "Run" in 'Editor' tab.
% Or section per section by clicking on "Run section"

%% Input data
%input filename(s)
file = 'VP GABA - grp6';
experiment = 'strawberry milk drinking';
folder = 'Z:\Hnasko lab_Lauren Faget_Data\Hnasko lab_Fiber Photometry\RESUME_FPData_VP GCAMP_cohort1&2&3';

%input notes for saved data
note_condition1 = 'all licks';
note_condition2 = 'licks 3s ILI';
BOUT_TIME_THRESHOLD = 3; % example bout time threshold, in seconds

%% MODIFY trigger names here
%trigger names can be found in data/ epocs/ 
%format is '3 letters - underscore' or 'x - 3 numbers/numbers-letters combo - underscore'

EPOC_1 = 'RRli'; % enter name of trigger 1
EPOC_2 = 'LLli'; % enter name of trigger 1
EPOC_3 = 'Rbo_'; % enter name of trigger 1
EPOC_4 = 'LLbo'; % enter name of trigger 1

%% Import data - Enter tank file paths here

fp1='C:\Users\laure\Documents\--Labos\Hnasko lab_Fiber photometry\RESUME_FPData_VP GCAMP_cohort1&2&3\Milkshake drinking assay\Strawberry milk drinking_data\Hnasko_lab-190604\B0014_VGAT-190604-113850';
data1 = TDTbin2mat(fp1);

fp2='C:\Users\laure\Documents\--Labos\Hnasko lab_Fiber Photometry\RESUME_FPData_VP GCAMP_cohort1&2&3\Milkshake drinking assay\Strawberry milk drinking_data\Hnasko_lab-190604\B0015_VGAT-190604-135043';
data2 = TDTbin2mat(fp2);

fp3='C:\Users\laure\Documents\--Labos\Hnasko lab_Fiber Photometry\RESUME_FPData_VP GCAMP_cohort1&2&3\Milkshake drinking assay\Strawberry milk drinking_data\Hnasko_lab-190604\B0017_VGAT-190604-142503';
data3 = TDTbin2mat(fp3);

fp4='C:\Users\laure\Documents\--Labos\Hnasko lab_Fiber Photometry\RESUME_FPData_VP GCAMP_cohort1&2&3\Milkshake drinking assay\Strawberry milk drinking_data\Hnasko_lab-191003\B1104_VGAT-191003-110332';
data4 = TDTbin2mat(fp4);

fp5='C:\Users\laure\Documents\--Labos\Hnasko lab_Fiber Photometry\RESUME_FPData_VP GCAMP_cohort1&2&3\Milkshake drinking assay\Strawberry milk drinking_data\Hnasko_lab-191003\B1105_VGAT-191003-115553';
data5 = TDTbin2mat(fp5);

fp6='C:\Users\laure\Documents\--Labos\Hnasko lab_Fiber Photometry\RESUME_FPData_VP GCAMP_cohort1&2&3\Milkshake drinking assay\Strawberry milk drinking_data\Hnasko_lab-191003\B1106_VGAT-191003-140221';
data6 = TDTbin2mat(fp6);

%% 
% To call the .mat file instead of the data tank:
% load 'filepath.mat';

% See end of script for how to save .mat file

%% Calculate dF/F

%Gc = 465nm signal. initially named after GCAMP
%af = 405nm signal. autofluorescence
%cf = controlfit of 405. 405 aligned on 465 signal
% eeggilt function applies a low (x,x,low,x) and high (x,x,x,high) pass filter to data
%fs = TDT processor sampling rate - 1017.25 data point per second (1017.25Hz)

Gc1=eegfilt(data1.streams.x465F.data,data1.streams.x465F.fs,0,10); 
af1=eegfilt(data1.streams.x405F.data,data1.streams.x405F.fs,0,10); 
cf1=controlFit3(Gc1,af1);

Gc2=eegfilt(data2.streams.x465F.data,data2.streams.x465F.fs,0,10); 
af2=eegfilt(data2.streams.x405F.data,data2.streams.x405F.fs,0,10); 
cf2=controlFit3(Gc2,af2);

Gc3=eegfilt(data3.streams.x465F.data,data3.streams.x465F.fs,0,10); 
af3=eegfilt(data3.streams.x405F.data,data3.streams.x405F.fs,0,10); 
cf3=controlFit3(Gc3,af3);

Gc4=eegfilt(data4.streams.x465F.data,data4.streams.x465F.fs,0,10); 
af4=eegfilt(data4.streams.x405F.data,data4.streams.x405F.fs,0,10); 
cf4=controlFit3(Gc4,af4);

Gc5=eegfilt(data5.streams.x465F.data,data5.streams.x465F.fs,0,10); 
af5=eegfilt(data5.streams.x405F.data,data5.streams.x405F.fs,0,10); 
cf5=controlFit3(Gc5,af5);

Gc6=eegfilt(data6.streams.x465F.data,data6.streams.x465F.fs,0,10); 
af6=eegfilt(data6.streams.x405F.data,data6.streams.x405F.fs,0,10); 
cf6=controlFit3(Gc6,af6);

fs=data1.streams.x465F.fs;
 
[dF1] = deltaFF (Gc1, cf1);
[dF2] = deltaFF (Gc2, cf2);
[dF3] = deltaFF (Gc3, cf3);
[dF4] = deltaFF (Gc4, cf4);
[dF5] = deltaFF (Gc5, cf5);
[dF6] = deltaFF (Gc6, cf6);

%%  this will plot the non-normalized Gcamp/AF signals together
% figure;
% hold on;
% plot(Gc1,'k');
% plot(af1,'b');
% plot(cf1,'r');

%% B. Trigger timepoint extraction 
% plot the entire dF and the triggers you made

tr1_1=data1.epocs.(EPOC_1).onset;
tr1_2=data2.epocs.(EPOC_1).onset;
tr1_3=data3.epocs.(EPOC_1).onset;
tr1_4=data4.epocs.(EPOC_1).onset;
tr1_5=data5.epocs.(EPOC_1).onset;
tr1_6=data6.epocs.(EPOC_1).onset;

%
tr2_1=data1.epocs.(EPOC_3).onset;
tr2_2=data2.epocs.(EPOC_3).onset;
tr2_3=data3.epocs.(EPOC_3).onset;
tr2_4=data4.epocs.(EPOC_3).onset;
tr2_5=data5.epocs.(EPOC_3).onset;
tr2_6=data6.epocs.(EPOC_3).onset;

% 
trt1_1=(tr1_1*fs);
trt1_2=(tr1_2*fs);
trt1_3=(tr1_3*fs);
trt1_4=(tr1_4*fs);
trt1_5=(tr1_5*fs);
trt1_6=(tr1_6*fs);

%
trt2_1=(tr2_1*fs);
trt2_2=(tr2_2*fs);
trt2_3=(tr2_3*fs);
trt2_4=(tr2_4*fs);
trt2_5=(tr2_5*fs);
trt2_6=(tr2_6*fs);

%% %% %% Create a matrix of LICK onset depending on ILI 


lick_on_diffA = diff(tr1_1);
lick_diff_indA = find(lick_on_diffA >= BOUT_TIME_THRESHOLD);
% 
diff_ind_iA = 1;
for iA = 1:length(lick_diff_indA)
    data1.epocs.Lbout.onset(iA) = tr1_1(diff_ind_iA);
    data1.epocs.Lbout.data(iA) = 1; % set the data value, arbitrary 1
    diff_ind_iA = lick_diff_indA(iA) + 1; % increment the index
end

lick_on_diffB = diff(tr1_2);
lick_diff_indB = find(lick_on_diffB >= BOUT_TIME_THRESHOLD);
% 
diff_ind_iB = 1;
for iB = 1:length(lick_diff_indB)
    data2.epocs.Lbout.onset(iB) = tr1_2(diff_ind_iB);
    data2.epocs.Lbout.data(iB) = 1; % set the data value, arbitrary 1
    diff_ind_iB = lick_diff_indB(iB) + 1; % increment the index
end

lick_on_diffC = diff(tr1_3);
lick_diff_indC = find(lick_on_diffC >= BOUT_TIME_THRESHOLD);
% 
diff_ind_iC = 1;
for iC = 1:length(lick_diff_indC)
    data3.epocs.Lbout.onset(iC) = tr1_3(diff_ind_iC);
    data3.epocs.Lbout.data(iC) = 1; % set the data value, arbitrary 1
    diff_ind_iC = lick_diff_indC(iC) + 1; % increment the index
end

lick_on_diffD = diff(tr1_4);
lick_diff_indD = find(lick_on_diffD >= BOUT_TIME_THRESHOLD);
% 
diff_ind_iD = 1;
for iD = 1:length(lick_diff_indD)
    data4.epocs.Lbout.onset(iD) = tr1_4(diff_ind_iD);
    data4.epocs.Lbout.data(iD) = 1; % set the data value, arbitrary 1
    diff_ind_iD = lick_diff_indD(iD) + 1; % increment the index
end

lick_on_diffE = diff(tr1_5);
lick_diff_indE = find(lick_on_diffE >= BOUT_TIME_THRESHOLD);
% 
diff_ind_iE = 1;
for iE = 1:length(lick_diff_indE)
    data5.epocs.Lbout.onset(iE) = tr1_5(diff_ind_iE);
    data5.epocs.Lbout.data(iE) = 1; % set the data value, arbitrary 1
    diff_ind_iE = lick_diff_indE(iE) + 1; % increment the index
end

lick_on_diffF = diff(tr1_6);
lick_diff_indF = find(lick_on_diffF >= BOUT_TIME_THRESHOLD);
% 
diff_ind_iF = 1;
for iF = 1:length(lick_diff_indF)
    data6.epocs.Lbout.onset(iF) = tr1_6(diff_ind_iF);
    data6.epocs.Lbout.data(iF) = 1; % set the data value, arbitrary 1
    diff_ind_iF = lick_diff_indF(iF) + 1; % increment the index
end


% Transpose the arrays to make them column vectors like other epocs
data1.epocs.Lbout.onset = tr1_1([1; lick_diff_indA+1])';
data2.epocs.Lbout.onset = tr1_2([1; lick_diff_indB+1])';
data3.epocs.Lbout.onset = tr1_3([1; lick_diff_indC+1])';
data4.epocs.Lbout.onset = tr1_4([1; lick_diff_indD+1])';
data5.epocs.Lbout.onset = tr1_5([1; lick_diff_indE+1])';
data6.epocs.Lbout.onset = tr1_6([1; lick_diff_indF+1])';

%% Create new trigger for ILI-dependent lick onsets

tr3_1=data1.epocs.Lbout.onset';
tr3_2=data2.epocs.Lbout.onset';
tr3_3=data3.epocs.Lbout.onset';
tr3_4=data4.epocs.Lbout.onset';
tr3_5=data5.epocs.Lbout.onset';
tr3_6=data6.epocs.Lbout.onset';

trt3_1=(tr3_1*fs);
trt3_2=(tr3_2*fs);
trt3_3=(tr3_3*fs);
trt3_4=(tr3_4*fs);
trt3_5=(tr3_5*fs);
trt3_6=(tr3_6*fs);


%% Plot all lick onsets, medpc bout onsets, and ILI-dependent lick onsets

figure;
title('Entire Experiment - animal 1');hold on;
plot((1/fs:1/fs:length(dF1)/fs),dF1,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_1)
     plot([tr1_1(i);tr1_1(i)],[-15;15],'c');
 end
  for i = 1:length(tr2_1)
     plot([tr2_1(i);tr2_1(i)],[-15;15],'b');
 end
  for i = 1:length(tr3_1)
     plot([tr3_1(i);tr3_1(i)],[-15;15],'r');
 end

 figure;
title('Entire Experiment - animal 2');hold on;
plot((1/fs:1/fs:length(dF2)/fs),dF2,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

 for i = 1:length(tr1_2)
     plot([tr1_2(i);tr1_2(i)],[-15;15],'c');
 end
  for i = 1:length(tr2_2)
     plot([tr2_2(i);tr2_2(i)],[-15;15],'b');
 end
  for i = 1:length(tr3_2)
     plot([tr3_2(i);tr3_2(i)],[-15;15],'r');
 end

figure;
title('Entire Experiment - animal 3');hold on;
plot((1/fs:1/fs:length(dF3)/fs),dF3,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_3)
     plot([tr1_3(i);tr1_3(i)],[-15;15],'c');
 end
   for i = 1:length(tr2_3)
     plot([tr2_3(i);tr2_3(i)],[-15;15],'b');
 end
  for i = 1:length(tr3_3)
     plot([tr3_3(i);tr3_3(i)],[-15;15],'r');
 end

figure;
title('Entire Experiment - animal 4');hold on;
plot((1/fs:1/fs:length(dF4)/fs),dF4,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_4)
     plot([tr1_4(i);tr1_4(i)],[-15;15],'c');
 end
  for i = 1:length(tr2_4)
     plot([tr2_4(i);tr2_4(i)],[-15;15],'b');
 end
  for i = 1:length(tr3_4)
     plot([tr3_4(i);tr3_4(i)],[-15;15],'r');
 end

figure;
title('Entire Experiment - animal 5');hold on;
plot((1/fs:1/fs:length(dF5)/fs),dF5,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_5)
     plot([tr1_5(i);tr1_5(i)],[-15;15],'c');
 end
  for i = 1:length(tr2_5)
     plot([tr2_5(i);tr2_5(i)],[-15;15],'b');
 end
  for i = 1:length(tr3_5)
     plot([tr3_5(i);tr3_5(i)],[-15;15],'r');
 end

 figure;
title('Entire Experiment - animal 6');hold on;
plot((1/fs:1/fs:length(dF6)/fs),dF6,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_6)
     plot([tr1_6(i);tr1_6(i)],[-15;15],'c');
 end
  for i = 1:length(tr2_6)
     plot([tr2_6(i);tr2_6(i)],[-15;15],'b');
 end
  for i = 1:length(tr3_6)
     plot([tr3_6(i);tr3_6(i)],[-15;15],'r');
 end

% close all;

%% C.this plots +/- time windows around the triggers made during experiment.

%DEFINE your time window here;
bltimewin=5;  %% baseline time window in seconds
timewin=10;   %% time window after trigger in seconds
timescale = bltimewin + timewin;

dF_trials1_1=[]; 
dF_trials1_2=[]; 
dF_trials1_3=[]; 
dF_trials1_4=[]; 
dF_trials1_5=[]; 
dF_trials1_6=[]; 

dF_trials3_1=[]; 
dF_trials3_2=[]; 
dF_trials3_3=[]; 
dF_trials3_4=[]; 
dF_trials3_5=[]; 
dF_trials3_6=[]; 


%animal 1
% extract trials for the first trigger.
count=0;
for i = 1:length(trt1_1)
     if (trt1_1(i)+timewin*fs<length(dF1) && trt1_1(i)-timewin*fs>0)     
        count=count+1;
        dF_trials1_1(count,:)=dF1(trt1_1(i)-(bltimewin*fs):trt1_1(i)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i = 1:length(trt3_1)
     if (trt3_1(i)+timewin*fs<length(dF1) && trt3_1(i)-timewin*fs>0)     
        count=count+1;
        dF_trials3_1(count,:)=dF1(trt3_1(i)-(bltimewin*fs):trt3_1(i)+(timewin*fs));  
    end
end

%animal 2
% extract trials for the first trigger.
count=0;
for i = 1:length(trt1_2)
     if (trt1_2(i)+timewin*fs<length(dF2) && trt1_2(i)-timewin*fs>0)     
        count=count+1;
        dF_trials1_2(count,:)=dF2(trt1_2(i)-(bltimewin*fs):trt1_2(i)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i = 1:length(trt3_2)
     if (trt3_2(i)+timewin*fs<length(dF2) && trt3_2(i)-timewin*fs>0)     
        count=count+1;
        dF_trials3_2(count,:)=dF2(trt3_2(i)-(bltimewin*fs):trt3_2(i)+(timewin*fs));  
    end
end

%animal 3
% extract trials for the first trigger.
count=0;
for i = 1:length(trt1_3)
     if (trt1_3(i)+timewin*fs<length(dF3) && trt1_3(i)-timewin*fs>0)     
        count=count+1;
        dF_trials1_3(count,:)=dF3(trt1_3(i)-(bltimewin*fs):trt1_3(i)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i = 1:length(trt3_3)
     if (trt3_3(i)+timewin*fs<length(dF3) && trt3_3(i)-timewin*fs>0)     
        count=count+1;
        dF_trials3_3(count,:)=dF3(trt3_3(i)-(bltimewin*fs):trt3_3(i)+(timewin*fs));  
    end
end

%animal 4
% extract trials for the first trigger.
count=0;
for i = 1:length(trt1_4)
     if (trt1_4(i)+timewin*fs<length(dF4) && trt1_4(i)-timewin*fs>0)     
        count=count+1;
        dF_trials1_4(count,:)=dF4(trt1_4(i)-(bltimewin*fs):trt1_4(i)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i = 1:length(trt3_4)
     if (trt3_4(i)+timewin*fs<length(dF4) && trt3_4(i)-timewin*fs>0)     
        count=count+1;
        dF_trials3_4(count,:)=dF4(trt3_4(i)-(bltimewin*fs):trt3_4(i)+(timewin*fs));  
    end
end

%animal 5
% extract trials for the first trigger.
count=0;
for i = 1:length(trt1_5)
     if (trt1_5(i)+timewin*fs<length(dF5) && trt1_5(i)-timewin*fs>0)     
        count=count+1;
        dF_trials1_5(count,:)=dF5(trt1_5(i)-(bltimewin*fs):trt1_5(i)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i = 1:length(trt3_5)
     if (trt3_5(i)+timewin*fs<length(dF5) && trt3_5(i)-timewin*fs>0)     
        count=count+1;
        dF_trials3_5(count,:)=dF5(trt3_5(i)-(bltimewin*fs):trt3_5(i)+(timewin*fs));  
    end
end

%animal 6
% extract trials for the third trigger.
count=0;
for i = 1:length(trt1_6)
     if (trt1_6(i)+timewin*fs<length(dF6) && trt1_6(i)-timewin*fs>0)     
        count=count+1;
        dF_trials1_6(count,:)=dF6(trt1_6(i)-(bltimewin*fs):trt1_6(i)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i = 1:length(trt3_6)
     if (trt3_6(i)+timewin*fs<length(dF6) && trt3_6(i)-timewin*fs>0)     
        count=count+1;
        dF_trials3_6(count,:)=dF6(trt3_6(i)-(bltimewin*fs):trt3_6(i)+(timewin*fs));  
    end
end


%Non-normalized dF/F for all trials and all animals
dF_trials1= [dF_trials1_1;dF_trials1_2;dF_trials1_3;dF_trials1_4;dF_trials1_5;dF_trials1_6];
dF_trials3= [dF_trials3_1;dF_trials3_2;dF_trials3_3;dF_trials3_4;dF_trials3_5;dF_trials3_6];


%%  extract Normalized data per trigger and graph ColorBars

% %VALUES TO MODIFY ACCORDING TO NEED HERE
% A=(bltimewin-1)*fs; % beginning of baseline to normalize with
% B=bltimewin*fs; % end of baseline to normalize with
% % time unit is in ms. Multiply value by fs to obtain seconds.
% % normalize data using a baseline 1 second prior to trigger in this
% % example
% 
% BASELINE_NORM = [A B];
% 
% bl1_1=dF_trials1_1(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
% blm1_1=mean(bl1_1,2);
% dFn1_1=(dF_trials1_1-repmat(blm1_1,1,size(dF_trials1_1,2)));
% 
% bl1_2=dF_trials1_2(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
% blm1_2=mean(bl1_2,2);
% dFn1_2=(dF_trials1_2-repmat(blm1_2,1,size(dF_trials1_2,2)));
% 
% bl1_3=dF_trials1_3(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
% blm1_3=mean(bl1_3,2);
% dFn1_3=(dF_trials1_3-repmat(blm1_3,1,size(dF_trials1_3,2)));
% 
% bl1_4=dF_trials1_4(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
% blm1_4=mean(bl1_4,2);
% dFn1_4=(dF_trials1_4-repmat(blm1_4,1,size(dF_trials1_4,2)));
% 
% bl1_5=dF_trials1_5(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
% blm1_5=mean(bl1_5,2);
% dFn1_5=(dF_trials1_5-repmat(blm1_5,1,size(dF_trials1_5,2)));
% 
% bl1_6=dF_trials1_6(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
% blm1_6=mean(bl1_6,2);
% dFn1_6=(dF_trials1_6-repmat(blm1_6,1,size(dF_trials1_6,2)));
% 
% 
% %Normalized dF/F for all trials and all animals
% dFn1= [dFn1_1;dFn1_2;dFn1_3;dFn1_4;dFn1_5;dFn1_6];
% 

% figure;
% imagesc(dFn1); colorbar; hold on;
% caxis([-4 8]);  % Y axis scale
% xline(5*fs,'-k','linewidth',2); % create a line at trigger onset
% xticks((0:fs:timescale*fs));  % X axis scale
% xticklabels({-5 : 10});  % X axis scale labeling
% title('normalized dF/F trigger 1');
% ylabel('Trials', 'FontSize', 12);
% xlabel('Time (s)','FontSize',12)
% set (gca,'TickDir','out');
% 
% 
%% Plot several graphs in one figure
% figure;
% subplot(221);
% imagesc(dFn1); colorbar; title('Normalized data trigger 1');
% caxis([-5 20]);
% xticks((0:fs:timescale*fs));
% xticklabels({-5 : 10}); hold on;
% subplot(222);
% imagesc(dF_trials1);colorbar;title('Non-Normalized data trigger 1');
% caxis([-5 20]);
% xticks((0:fs:timescale*fs));
% xticklabels({-5 : 10}); hold on;


%% Normalized Mean +/- SEM with shadded error boars

% dFmnnm1=mean(dFn1(1:end,1:end));
% dFsemnm1=std(dFn1(1:end,1:end))/sqrt(size(dFn1,1));
% 
% figure;
% shadedErrorBar((-bltimewin:1/fs:timewin),dFmnnm1,dFsemnm1); 
% title('normalized dF/F trigger 1');
% axis([-bltimewin inf -4 8]);
% hold on;
% plot([0;0],[-4;8],'b','Linewidth',1);
% plot([5;5],[-4;8],':b','Linewidth',1); % plot a second line at the end of the stim
% ylabel('dF/F', 'FontSize', 12);
% xlabel('Time (s)','FontSize',12);
% set (gca,'TickDir','out');

 %% plotting triggers against each other normalized
% figure;
% hold on;
% plot((-bltimewin:1/fs:timewin),dFmnnm1,'c');
% plot((-bltimewin:1/fs:timewin),dFmnnm2,'g');
% plot((-bltimewin:1/fs:timewin),dFmnnm3,'r');
% plot((-bltimewin:1/fs:timewin),dFmnnm4,'b');
% axis([-bltimewin inf -5 20]);
% hold on;
% plot([0;0],[-5;20],'k','Linewidth',1);
% plot([1;1],[-5;20],':k','Linewidth',1);
% legend('trigger 1','trigger 2','trigger 3','trigger 4');
% title('Normalized DF/F laser');
% ylabel('dF/F', 'FontSize', 12);
% xlabel('Time (s)','FontSize',12);
% set (gca,'TickDir','out');

%%  Merge days
% 
% dF_trials1 = [dF_trials1a; dF_trials1b];
% 
% dF_trials1_1 = [dF_trials1_1a; dF_trials1_1b];
% dF_trials1_2 = [dF_trials1_2a; dF_trials1_2b];
% dF_trials1_4 = [dF_trials1_4a; dF_trials1_4b];
% dF_trials1_6 = [dF_trials1_6a; dF_trials1_6b];
% 
% 
%% Z SCORE at lick and ILI-dependent lick onset

%VALUES TO MODIFY ACCORDING TO NEED HERE
C=1; % beginning of baseline  (1 means first data point of the 5s baseline)
D=bltimewin*fs; % end of baseline
% time unit is in ms. Multiply value by fs to obtain seconds.
% This example use a baseline of 5 seconds prior to trigger

BASELINE_Z = [C D]; 

%TRIGGER 1
zall1 = zeros(size(dF_trials1));
for i = 1:size(dF_trials1,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials1(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials1(i,ind)); % baseline period stdev
    zall1(i,:)=(dF_trials1(i,:) - zb)/zsd; % Z score per bin
end
zerror1 = std(zall1)/sqrt(size(zall1,1));
zmean1 =  mean (zall1);

%FIGURES
% heat map - colorbar graph - trigger 3
figure; 
imagesc(zall1); colorbar; hold on;
caxis([-4 6]);
xticks((0:fs:timescale*fs));
xticklabels({-5 : 10});
title('Z-Score Trigger 1 Heat Map')
ylabel('Trials', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition1,'-','zscore colorplot'); 
title(note_title);

% save figure as matlab fig
% cd (folder);
% saveas(gcf,note_title, 'fig');

% mean +/- sem shadded error bar graph - trigger 1
figure; 
shadedErrorBar([-bltimewin:1/fs:timewin],zmean1,zerror1, 'lineProps','b'); hold on;
axis([-5 inf -4 6]);
plot([0;0],[-4;6],'k','Linewidth',1);
ylabel('z score', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition1,'-','zscore'); 
title(note_title);

% save figure as matlab fig
% cd (folder);
% saveas(gcf,note_title, 'fig');


%TRIGGER 3
zall3 = zeros(size(dF_trials3));
for i = 1:size(dF_trials3,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3(i,ind)); % baseline period stdev
    zall3(i,:)=(dF_trials3(i,:) - zb)/zsd; % Z score per bin
end
zerror3 = std(zall3)/sqrt(size(zall3,1));
zmean3 =  mean (zall3);

%FIGURES
% heat map - colorbar graph - trigger 3
figure; 
imagesc(zall3); colorbar; hold on;
caxis([-4 6]);
xticks((0:fs:timescale*fs));
xticklabels({-5 : 10});
title('Z-Score Trigger 1 Heat Map')
ylabel('Trials', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition2,'-','zscore colorplot'); 
title(note_title);

% save figure as matlab fig
% cd (folder);
% saveas(gcf,note_title, 'fig');

% mean +/- sem shadded error bar graph - trigger 3
figure; 
shadedErrorBar([-bltimewin:1/fs:timewin],zmean3,zerror3, 'lineProps','b'); hold on;
axis([-5 inf -4 6]);
plot([0;0],[-4;6],'k','Linewidth',1);
ylabel('z score', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition2,'-','zscore'); 
title(note_title);

% save figure as matlab fig
% cd (folder);
% saveas(gcf,note_title, 'fig');
% 
%%  Quantify min/max (peak) z score value for ILI-dependent lick onset
% per 1s bin

E=(bltimewin-1)*fs; 
F=bltimewin*fs; 
G=(bltimewin+1)*fs;
H=(bltimewin+4)*fs;
I=(bltimewin+5)*fs;
J=(bltimewin+6)*fs;
K=(bltimewin+2)*fs;
L=(bltimewin-2)*fs;
M=(bltimewin+7)*fs;
N=(bltimewin+3)*fs;

%TRIGGER 3  - animal 1
zall3_1 = zeros(size(dF_trials3_1));
for i = 1:size(dF_trials3_1,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_1(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_1(i,ind)); % baseline period stdev
    zall3_1(i,:)=(dF_trials3_1(i,:) - zb)/zsd; % Z score per bin
end
zerror3_1 = std(zall3_1)/sqrt(size(zall3_1,1));
zmean3_1a =  mean (zall3_1);
zmean3_1 = movmean(zmean3_1a, 10);

zmean3_1A= zmean3_1 (:,E:F);
zmean3_1B= zmean3_1 (:,F:G);
zmean3_1C= zmean3_1 (:,H:I);
zmean3_1D= zmean3_1 (:,I:J);
zmean3_1E= zmean3_1 (:,G:K);
zmean3_1F= zmean3_1 (:,L:E);
zmean3_1G= zmean3_1 (:,J:M);
zmean3_1H= zmean3_1 (:,K:N);
zmean3_1I= zmean3_1 (:,N:H);

[~,X] = max(abs(zmean3_1A));
zmean3_1Apeak = zmean3_1A(X);
[~,X] = max(abs(zmean3_1B));
zmean3_1Bpeak = zmean3_1B(X);
[~,X] = max(abs(zmean3_1C));
zmean3_1Cpeak = zmean3_1C(X);
[~,X] = max(abs(zmean3_1D));
zmean3_1Dpeak = zmean3_1D(X);
[~,X] = max(abs(zmean3_1E));
zmean3_1Epeak = zmean3_1E(X);
[~,X] = max(abs(zmean3_1F));
zmean3_1Fpeak = zmean3_1F(X);
[~,X] = max(abs(zmean3_1G));
zmean3_1Gpeak = zmean3_1G(X);
[~,X] = max(abs(zmean3_1H));
zmean3_1Hpeak = zmean3_1H(X);
[~,X] = max(abs(zmean3_1I));
zmean3_1Ipeak = zmean3_1I(X);


%TRIGGER 3  - animal 2
zall3_2 = zeros(size(dF_trials3_2));
for i = 1:size(dF_trials3_2,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_2(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_2(i,ind)); % baseline period stdev
    zall3_2(i,:)=(dF_trials3_2(i,:) - zb)/zsd; % Z score per bin
end
zerror3_2 = std(zall3_2)/sqrt(size(zall3_2,1));
zmean3_2a =  mean (zall3_2);
zmean3_2 = movmean(zmean3_2a, 10);

zmean3_2A= zmean3_2 (:,E:F);
zmean3_2B= zmean3_2 (:,F:G);
zmean3_2C= zmean3_2 (:,H:I);
zmean3_2D= zmean3_2 (:,I:J);
zmean3_2E= zmean3_2 (:,G:K);
zmean3_2F= zmean3_2 (:,L:E);
zmean3_2G= zmean3_2 (:,J:M);
zmean3_2H= zmean3_2 (:,K:N);
zmean3_2I= zmean3_2 (:,N:H);

[~,X] = max(abs(zmean3_2A));
zmean3_2Apeak = zmean3_2A(X);
[~,X] = max(abs(zmean3_2B));
zmean3_2Bpeak = zmean3_2B(X);
[~,X] = max(abs(zmean3_2C));
zmean3_2Cpeak = zmean3_2C(X);
[~,X] = max(abs(zmean3_2D));
zmean3_2Dpeak = zmean3_2D(X);
[~,X] = max(abs(zmean3_2E));
zmean3_2Epeak = zmean3_2E(X);
[~,X] = max(abs(zmean3_2F));
zmean3_2Fpeak = zmean3_2F(X);
[~,X] = max(abs(zmean3_2G));
zmean3_2Gpeak = zmean3_2G(X);
[~,X] = max(abs(zmean3_2H));
zmean3_2Hpeak = zmean3_2H(X);
[~,X] = max(abs(zmean3_2I));
zmean3_2Ipeak = zmean3_2I(X);


%TRIGGER 3  - animal 3
zall3_3 = zeros(size(dF_trials3_3));
for i = 1:size(dF_trials3_3,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_3(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_3(i,ind)); % baseline period stdev
    zall3_3(i,:)=(dF_trials3_3(i,:) - zb)/zsd; % Z score per bin
end
zerror3_3 = std(zall3_3)/sqrt(size(zall3_3,1));
zmean3_3a =  mean (zall3_3);
zmean3_3 = movmean(zmean3_3a, 10);

zmean3_3A= zmean3_3 (:,E:F);
zmean3_3B= zmean3_3 (:,F:G);
zmean3_3C= zmean3_3 (:,H:I);
zmean3_3D= zmean3_3 (:,I:J);
zmean3_3E= zmean3_3 (:,G:K);
zmean3_3F= zmean3_3 (:,L:E);
zmean3_3G= zmean3_3 (:,J:M);
zmean3_3H= zmean3_3 (:,K:N);
zmean3_3I= zmean3_3 (:,N:H);

[~,X] = max(abs(zmean3_3A));
zmean3_3Apeak = zmean3_3A(X);
[~,X] = max(abs(zmean3_3B));
zmean3_3Bpeak = zmean3_3B(X);
[~,X] = max(abs(zmean3_3C));
zmean3_3Cpeak = zmean3_3C(X);
[~,X] = max(abs(zmean3_3D));
zmean3_3Dpeak = zmean3_3D(X);
[~,X] = max(abs(zmean3_3E));
zmean3_3Epeak = zmean3_3E(X);
[~,X] = max(abs(zmean3_3F));
zmean3_3Fpeak = zmean3_3F(X);
[~,X] = max(abs(zmean3_3G));
zmean3_3Gpeak = zmean3_3G(X);
[~,X] = max(abs(zmean3_3H));
zmean3_3Hpeak = zmean3_3H(X);
[~,X] = max(abs(zmean3_3I));
zmean3_3Ipeak = zmean3_3I(X);


%TRIGGER 3  - animal 4
zall3_4 = zeros(size(dF_trials3_4));
for i = 1:size(dF_trials3_4,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_4(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_4(i,ind)); % baseline period stdev
    zall3_4(i,:)=(dF_trials3_4(i,:) - zb)/zsd; % Z score per bin
end
zerror3_4 = std(zall3_4)/sqrt(size(zall3_4,1));
zmean3_4a =  mean (zall3_4);
zmean3_4 = movmean(zmean3_4a, 10);

zmean3_4A= zmean3_4 (:,E:F);
zmean3_4B= zmean3_4 (:,F:G);
zmean3_4C= zmean3_4 (:,H:I);
zmean3_4D= zmean3_4 (:,I:J);
zmean3_4E= zmean3_4 (:,G:K);
zmean3_4F= zmean3_4 (:,L:E);
zmean3_4G= zmean3_4 (:,J:M);
zmean3_4H= zmean3_4 (:,K:N);
zmean3_4I= zmean3_4 (:,N:H);

[~,X] = max(abs(zmean3_4A));
zmean3_4Apeak = zmean3_4A(X);
[~,X] = max(abs(zmean3_4B));
zmean3_4Bpeak = zmean3_4B(X);
[~,X] = max(abs(zmean3_4C));
zmean3_4Cpeak = zmean3_4C(X);
[~,X] = max(abs(zmean3_4D));
zmean3_4Dpeak = zmean3_4D(X);
[~,X] = max(abs(zmean3_4E));
zmean3_4Epeak = zmean3_4E(X);
[~,X] = max(abs(zmean3_4F));
zmean3_4Fpeak = zmean3_4F(X);
[~,X] = max(abs(zmean3_4G));
zmean3_4Gpeak = zmean3_4G(X);
[~,X] = max(abs(zmean3_4H));
zmean3_4Hpeak = zmean3_4H(X);
[~,X] = max(abs(zmean3_4I));
zmean3_4Ipeak = zmean3_4I(X);


% %TRIGGER 3  - animal 5
zall3_5 = zeros(size(dF_trials3_5));
for i = 1:size(dF_trials3_5,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_5(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_5(i,ind)); % baseline period stdev
    zall3_5(i,:)=(dF_trials3_5(i,:) - zb)/zsd; % Z score per bin
end
zerror3_5 = std(zall3_5)/sqrt(size(zall3_5,1));
zmean3_5a =  mean (zall3_5);
zmean3_5 = movmean(zmean3_5a, 10);
% 
zmean3_5A= zmean3_5 (:,E:F);
zmean3_5B= zmean3_5 (:,F:G);
zmean3_5C= zmean3_5 (:,H:I);
zmean3_5D= zmean3_5 (:,I:J);
zmean3_5E= zmean3_5 (:,G:K);
zmean3_5F= zmean3_5 (:,L:E);
zmean3_5G= zmean3_5 (:,J:M);
zmean3_5H= zmean3_5 (:,K:N);
zmean3_5I= zmean3_5 (:,N:H);

 
[~,X] = max(abs(zmean3_5A));
zmean3_5Apeak = zmean3_5A(X);
[~,X] = max(abs(zmean3_5B));
zmean3_5Bpeak = zmean3_5B(X);
[~,X] = max(abs(zmean3_5C));
zmean3_5Cpeak = zmean3_5C(X);
[~,X] = max(abs(zmean3_5D));
zmean3_5Dpeak = zmean3_5D(X);
[~,X] = max(abs(zmean3_5E));
zmean3_5Epeak = zmean3_5E(X);
[~,X] = max(abs(zmean3_5F));
zmean3_5Fpeak = zmean3_5F(X);
[~,X] = max(abs(zmean3_5G));
zmean3_5Gpeak = zmean3_5G(X);
[~,X] = max(abs(zmean3_5H));
zmean3_5Hpeak = zmean3_5H(X);
[~,X] = max(abs(zmean3_5I));
zmean3_5Ipeak = zmean3_5I(X);

%
% %TRIGGER 3  - animal 6
zall3_6 = zeros(size(dF_trials3_6));
for i = 1:size(dF_trials3_6,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_6(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_6(i,ind)); % baseline period stdev
    zall3_6(i,:)=(dF_trials3_6(i,:) - zb)/zsd; % Z score per bin
end
zerror3_6 = std(zall3_6)/sqrt(size(zall3_6,1));
zmean3_6a =  mean (zall3_6);
zmean3_6 = movmean(zmean3_6a, 10);
% 
zmean3_6A= zmean3_6 (:,E:F);
zmean3_6B= zmean3_6 (:,F:G);
zmean3_6C= zmean3_6 (:,H:I);
zmean3_6D= zmean3_6 (:,I:J);
zmean3_6E= zmean3_6 (:,G:K);
zmean3_6F= zmean3_6 (:,L:E);
zmean3_6G= zmean3_6 (:,J:M);
zmean3_6H= zmean3_6 (:,K:N);
zmean3_6I= zmean3_6 (:,N:H);

[~,X] = max(abs(zmean3_6A));
zmean3_6Apeak = zmean3_6A(X);
[~,X] = max(abs(zmean3_6B));
zmean3_6Bpeak = zmean3_6B(X);
[~,X] = max(abs(zmean3_6C));
zmean3_6Cpeak = zmean3_6C(X);
[~,X] = max(abs(zmean3_6D));
zmean3_6Dpeak = zmean3_6D(X);
[~,X] = max(abs(zmean3_6E));
zmean3_6Epeak = zmean3_6E(X);
[~,X] = max(abs(zmean3_6F));
zmean3_6Fpeak = zmean3_6F(X);
[~,X] = max(abs(zmean3_6G));
zmean3_6Gpeak = zmean3_6G(X);
[~,X] = max(abs(zmean3_6H));
zmean3_6Hpeak = zmean3_6H(X);
[~,X] = max(abs(zmean3_6I));
zmean3_6Ipeak = zmean3_6I(X);

%
 
zabsmax_PreStim21= [zmean3_1Fpeak;zmean3_2Fpeak;zmean3_3Fpeak;zmean3_4Fpeak;zmean3_5Fpeak;zmean3_6Fpeak];
zabsmax_Prestim10= [zmean3_1Apeak;zmean3_2Apeak;zmean3_3Apeak;zmean3_4Apeak;zmean3_5Apeak;zmean3_6Apeak];
zabsmax_PostStim01= [zmean3_1Bpeak;zmean3_2Bpeak;zmean3_3Bpeak;zmean3_4Bpeak;zmean3_5Bpeak;zmean3_6Bpeak];
zabsmax_PostStim12= [zmean3_1Epeak;zmean3_2Epeak;zmean3_3Epeak;zmean3_4Epeak;zmean3_5Epeak;zmean3_6Epeak];
zabsmax_PostStim23= [zmean3_1Hpeak;zmean3_2Hpeak;zmean3_3Hpeak;zmean3_4Hpeak;zmean3_5Hpeak;zmean3_6Hpeak];
zabsmax_PostStim34= [zmean3_1Ipeak;zmean3_2Ipeak;zmean3_3Ipeak;zmean3_4Ipeak;zmean3_5Ipeak;zmean3_6Ipeak];
zabsmax_PostStim45= [zmean3_1Cpeak;zmean3_2Cpeak;zmean3_3Cpeak;zmean3_4Cpeak;zmean3_5Cpeak;zmean3_6Cpeak];
zabsmax_PostStim56= [zmean3_1Dpeak;zmean3_2Dpeak;zmean3_3Dpeak;zmean3_4Dpeak;zmean3_5Dpeak;zmean3_6Dpeak];
zabsmax_PostStim67= [zmean3_1Gpeak;zmean3_2Gpeak;zmean3_3Gpeak;zmean3_4Gpeak;zmean3_5Gpeak;zmean3_6Gpeak];
%
%%  Quantify changes as area under the curve - AUC

% %Modify time window as needed - example with 1s prior stim and 1s during
% 
% % TRIGGER 1
% AUC1=[]; % 
% AUC1(1,1)=trapz(zmean1(:,(bltimewin-1)*fs:(bltimewin)*fs));
% AUC1(1,2)=trapz(zmean1(:,(bltimewin)*fs:(bltimewin+1)*fs));
% 
% figure;
% hBar = bar(AUC1, 'FaceColor', [.8 .8 .8]);
% ylim([-500 2500]); hold on;
% % centers = get(hBar, 'XData');
% % plot(centers(1:2), [1 1]*AUC(1,2)*1.1, '-k', 'LineWidth', 2)
% % p1 = plot(mean(centers(1:2)), AUC(1,2)*1.2, '*k');
% set(gca,'xticklabel',{'baseline','stim'});
% set (gca,'TickDir','out');
% title({ 'AUC trigger 1'});
% 

%% % mean +/- sem shadded error bar graph - group
% figure; 
% shadedErrorBar([-bltimewin:1/fs:timewin],ILI1_zmean1,ILI1_zerror1, 'lineProps','c'); hold on;
% shadedErrorBar([-bltimewin:1/fs:timewin],ILI2_zmean1,ILI2_zerror1, 'lineProps','g'); hold on;
% shadedErrorBar([-bltimewin:1/fs:timewin],ILI5_zmean1,ILI5_zerror1, 'lineProps','b'); hold on;
% shadedErrorBar([-bltimewin:1/fs:timewin],ILI10_zmean1,ILI10_zerror1, 'lineProps','k'); hold on;
% shadedErrorBar([-bltimewin:1/fs:timewin],ILI20_zmean1,ILI20_zerror1, 'lineProps','r'); hold on;
% axis([-5 inf -2 4]);
% plot([0;0],[-2;4],'k','Linewidth',0.5);
% ylabel('z score', 'FontSize', 12);
% xlabel('Time (s)','FontSize',12);
% legend ('ILI 1s (cyan)','ILI 2s (green)','ILI 5s (blue)','ILI 10s (black)','ILI 20s (red)');
% set (gca,'TickDir','out');
% note_title =  strcat(file,'-',experiment,'-','all ILI','zscore'); 
% title(note_title);
% 
% % save figure as matlab fig
% % cd (folder);
% saveas(gcf,note_title, 'fig');


%% Save matlab data file
% You can save space by saving only necessary data - deleting unwanted data in the
% workspace - before saving .mat file
% 
% filename =  strcat(file,'-',experiment,'-',note_condition); 
% save (filename);
