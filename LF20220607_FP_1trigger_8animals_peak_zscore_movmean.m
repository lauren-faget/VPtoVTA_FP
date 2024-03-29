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
% Script can be run all at once by ckicking on "Run" in 'Editor' tab.
% Or section per section by clicking on "Run section"

%% Input data
%input filename(s)
file = 'VP GABA - VTA Glut-grp8';
experiment = 'ICSS 20s TO PPB';
folder = 'E:\Hnasko lab_Lauren Faget_Data\Hnasko lab_Fiber Photometry\RESUME_FP_VP ChR2-VTA GCAMP\z-score Fig 300dpi';

%input notes for saved data
note_condition = 'active nspk (no laser)';

%% MODIFY trigger names here
%trigger names can be found in data/ epocs/ 
%format is '3 letters - underscore' or 'x - 3 numbers/numbers-letters combo - underscore'

EPOC_1 = 'act_'; % enter name of trigger 1

%% Import data - Enter tank file paths here

fp1='E:\Hnasko lab_Lauren Faget_Data\Hnasko lab_Fiber Photometry\2021_Feb_VPgaba-ChR2_VTAglut-GCAMP6f_cohort2\ICSS\Hnasko_lab-210511\B5087_glut-210511-154450';
data1 = TDTbin2mat(fp1);

fp2='E:\Hnasko lab_Lauren Faget_Data\Hnasko lab_Fiber Photometry\2021_Feb_VPgaba-ChR2_VTAglut-GCAMP6f_cohort2\ICSS\Hnasko_lab-210511\B5089_glut-210511-145726';
data2 = TDTbin2mat(fp2);

fp3='E:\Hnasko lab_Lauren Faget_Data\Hnasko lab_Fiber Photometry\2021_Feb_VPgaba-ChR2_VTAglut-GCAMP6f_cohort2\ICSS\Hnasko_lab-210511\B5090_glut-210511-090957';
data3 = TDTbin2mat(fp3);

fp4='E:\Hnasko lab_Lauren Faget_Data\Hnasko lab_Fiber Photometry\2021_Feb_VPgaba-ChR2_VTAglut-GCAMP6f_cohort2\ICSS\Hnasko_lab-210511\B5091_glut-210511-104957';
data4 = TDTbin2mat(fp4);

fp5='E:\Hnasko lab_Lauren Faget_Data\Hnasko lab_Fiber Photometry\2021_Feb_VPgaba-ChR2_VTAglut-GCAMP6f_cohort2\ICSS\Hnasko_lab-210511\B5093_glut-210511-114324';
data5 = TDTbin2mat(fp5);

fp6='E:\Hnasko lab_Lauren Faget_Data\Hnasko lab_Fiber Photometry\2021_Feb_VPgaba-ChR2_VTAglut-GCAMP6f_cohort2\ICSS\Hnasko_lab-210511\B5096_glut-210511-123341';
data6 = TDTbin2mat(fp6);

fp7='E:\Hnasko lab_Lauren Faget_Data\Hnasko lab_Fiber Photometry\2021_Feb_VPgaba-ChR2_VTAglut-GCAMP6f_cohort2\ICSS\Hnasko_lab-210511\B5099_glut-210511-132148';
data7 = TDTbin2mat(fp7);

fp8='E:\Hnasko lab_Lauren Faget_Data\Hnasko lab_Fiber Photometry\2021_Feb_VPgaba-ChR2_VTAglut-GCAMP6f_cohort2\ICSS\Hnasko_lab-210511\B5100_glut-210511-140944';
data8 = TDTbin2mat(fp8);

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

Gc7=eegfilt(data7.streams.x465F.data,data7.streams.x465F.fs,0,10); 
af7=eegfilt(data7.streams.x405F.data,data7.streams.x405F.fs,0,10); 
cf7=controlFit3(Gc7,af7);

Gc8=eegfilt(data8.streams.x465F.data,data8.streams.x465F.fs,0,10); 
af8=eegfilt(data8.streams.x405F.data,data8.streams.x405F.fs,0,10); 
cf8=controlFit3(Gc8,af8);

fs=data1.streams.x465F.fs;
 
[dF1] = deltaFF (Gc1, cf1);
[dF2] = deltaFF (Gc2, cf2);
[dF3] = deltaFF (Gc3, cf3);
[dF4] = deltaFF (Gc4, cf4);
[dF5] = deltaFF (Gc5, cf5);
[dF6] = deltaFF (Gc6, cf6);
[dF7] = deltaFF (Gc7, cf7);
[dF8] = deltaFF (Gc8, cf8);

%% B. Trigger timepoint extraction 
% plot the entire dF and the triggers you made

tr1_1=data1.epocs.(EPOC_1).onset;

tr1_2=data2.epocs.(EPOC_1).onset;

tr1_3=data3.epocs.(EPOC_1).onset;

tr1_4=data4.epocs.(EPOC_1).onset;

tr1_5=data5.epocs.(EPOC_1).onset;

tr1_6=data6.epocs.(EPOC_1).onset;

tr1_7=data7.epocs.(EPOC_1).onset;

tr1_8=data8.epocs.(EPOC_1).onset;
% 

trt1_1=(tr1_1*fs);

trt1_2=(tr1_2*fs);

trt1_3=(tr1_3*fs);

trt1_4=(tr1_4*fs);

trt1_5=(tr1_5*fs);

trt1_6=(tr1_6*fs);

trt1_7=(tr1_7*fs);

trt1_8=(tr1_8*fs);


figure;
title('Entire Experiment - animal 1');hold on;
plot((1/fs:1/fs:length(dF1)/fs),dF1,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_1)
     plot([tr1_1(i);tr1_1(i)],[-15;15],'c');
 end

 figure;
title('Entire Experiment - animal 2');hold on;
plot((1/fs:1/fs:length(dF2)/fs),dF2,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

 for i = 1:length(tr1_2)
     plot([tr1_2(i);tr1_2(i)],[-15;15],'c');
 end

figure;
title('Entire Experiment - animal 3');hold on;
plot((1/fs:1/fs:length(dF3)/fs),dF3,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_3)
     plot([tr1_3(i);tr1_3(i)],[-15;15],'c');
 end
 
figure;
title('Entire Experiment - animal 4');hold on;
plot((1/fs:1/fs:length(dF4)/fs),dF4,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_4)
     plot([tr1_4(i);tr1_4(i)],[-15;15],'c');
 end

figure;
title('Entire Experiment - animal 5');hold on;
plot((1/fs:1/fs:length(dF5)/fs),dF5,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_5)
     plot([tr1_5(i);tr1_5(i)],[-15;15],'c');
 end

 figure;
title('Entire Experiment - animal 6');hold on;
plot((1/fs:1/fs:length(dF6)/fs),dF6,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_6)
     plot([tr1_6(i);tr1_6(i)],[-15;15],'c');
 end

  figure;
title('Entire Experiment - animal 7');hold on;
plot((1/fs:1/fs:length(dF7)/fs),dF7,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_7)
     plot([tr1_7(i);tr1_7(i)],[-15;15],'c');
 end

    figure;
title('Entire Experiment - animal 8');hold on;
plot((1/fs:1/fs:length(dF8)/fs),dF8,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_8)
     plot([tr1_8(i);tr1_8(i)],[-15;15],'c');
 end

close all;

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

dF_trials1_7=[]; 

dF_trials1_8=[]; 


%animal 1
% extract trials for the first trigger.
count=0;
for i = 1:length(trt1_1)
     if (trt1_1(i)+timewin*fs<length(dF1) && trt1_1(i)-timewin*fs>0)     
        count=count+1;
        dF_trials1_1(count,:)=dF1(trt1_1(i)-(bltimewin*fs):trt1_1(i)+(timewin*fs));  
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

%animal 3
% extract trials for the first trigger.
count=0;
for i = 1:length(trt1_3)
     if (trt1_3(i)+timewin*fs<length(dF3) && trt1_3(i)-timewin*fs>0)     
        count=count+1;
        dF_trials1_3(count,:)=dF3(trt1_3(i)-(bltimewin*fs):trt1_3(i)+(timewin*fs));  
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

%animal 5
% extract trials for the first trigger.
count=0;
for i = 1:length(trt1_5)
     if (trt1_5(i)+timewin*fs<length(dF5) && trt1_5(i)-timewin*fs>0)     
        count=count+1;
        dF_trials1_5(count,:)=dF5(trt1_5(i)-(bltimewin*fs):trt1_5(i)+(timewin*fs));  
    end
end

%animal 6
% extract trials for the first trigger.
count=0;
for i = 1:length(trt1_6)
     if (trt1_6(i)+timewin*fs<length(dF6) && trt1_6(i)-timewin*fs>0)     
        count=count+1;
        dF_trials1_6(count,:)=dF6(trt1_6(i)-(bltimewin*fs):trt1_6(i)+(timewin*fs));  
    end
end

%animal 7
% extract trials for the first trigger.
count=0;
for i = 1:length(trt1_7)
     if (trt1_7(i)+timewin*fs<length(dF7) && trt1_7(i)-timewin*fs>0)     
        count=count+1;
        dF_trials1_7(count,:)=dF7(trt1_7(i)-(bltimewin*fs):trt1_7(i)+(timewin*fs));  
    end
end

%animal 8
% extract trials for the first trigger.
count=0;
for i = 1:length(trt1_8)
     if (trt1_8(i)+timewin*fs<length(dF8) && trt1_8(i)-timewin*fs>0)     
        count=count+1;
        dF_trials1_8(count,:)=dF8(trt1_8(i)-(bltimewin*fs):trt1_8(i)+(timewin*fs));  
    end
end

%%  extract Normalized data per trigger and graph ColorBars

%VALUES TO MODIFY ACCORDING TO NEED HERE
A=(bltimewin-1)*fs; % beginning of baseline to normalize with
B=bltimewin*fs; % end of baseline to normalize with
% time unit is in ms. Multiply value by fs to obtain seconds.
% normalize data using a baseline 1 second prior to trigger in this
% example

BASELINE_NORM = [A B];

bl1_1=dF_trials1_1(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_1=mean(bl1_1,2);
dFn1_1=(dF_trials1_1-repmat(blm1_1,1,size(dF_trials1_1,2)));

bl1_2=dF_trials1_2(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_2=mean(bl1_2,2);
dFn1_2=(dF_trials1_2-repmat(blm1_2,1,size(dF_trials1_2,2)));

bl1_3=dF_trials1_3(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_3=mean(bl1_3,2);
dFn1_3=(dF_trials1_3-repmat(blm1_3,1,size(dF_trials1_3,2)));

bl1_4=dF_trials1_4(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_4=mean(bl1_4,2);
dFn1_4=(dF_trials1_4-repmat(blm1_4,1,size(dF_trials1_4,2)));

bl1_5=dF_trials1_5(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_5=mean(bl1_5,2);
dFn1_5=(dF_trials1_5-repmat(blm1_5,1,size(dF_trials1_5,2)));

bl1_6=dF_trials1_6(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_6=mean(bl1_6,2);
dFn1_6=(dF_trials1_6-repmat(blm1_6,1,size(dF_trials1_6,2)));

bl1_7=dF_trials1_7(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_7=mean(bl1_7,2);
dFn1_7=(dF_trials1_7-repmat(blm1_7,1,size(dF_trials1_7,2)));

bl1_8=dF_trials1_8(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_8=mean(bl1_8,2);
dFn1_8=(dF_trials1_8-repmat(blm1_8,1,size(dF_trials1_8,2)));


%Non-normalized dF/F for all trials and all animals
dF_trials1= [dF_trials1_1;dF_trials1_2;dF_trials1_3;dF_trials1_4;dF_trials1_5;dF_trials1_6;dF_trials1_7;dF_trials1_8];

%Normalized dF/F for all trials and all animals
dFn1= [dFn1_1;dFn1_2;dFn1_3;dFn1_4;dFn1_5;dFn1_6;dFn1_7;dFn1_8];


figure;
imagesc(dFn1); colorbar; hold on;
caxis([-12 20]);  % Y axis scale
xline(5*fs,'-k','linewidth',2); % create a line at trigger onset
xticks((0:fs:timescale*fs));  % X axis scale
xticklabels({-5 : 10});  % X axis scale labeling
title('normalized dF/F trigger 3');
ylabel('Trials', 'FontSize', 12);
xlabel('Time (s)','FontSize',12)
set (gca,'TickDir','out');


%% Normalized Mean +/- SEM with shadded error boars

dFmnnm1=mean(dFn1(1:end,1:end));
dFsemnm1=std(dFn1(1:end,1:end))/sqrt(size(dFn1,1));

figure;
shadedErrorBar((-bltimewin:1/fs:timewin),dFmnnm1,dFsemnm1); 
title('normalized dF/F trigger 3');
axis([-bltimewin inf -12 20]);
hold on;
plot([0;0],[-12;20],'b','Linewidth',1);
plot([5;5],[-12;20],':b','Linewidth',1); % plot a second line at the end of the stim
ylabel('dF/F', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');


%% Z SCORE

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
% heat map - colorbar graph - trigger 1
figure; 
imagesc(zall1); colorbar; hold on;
caxis([-12 20]);
xticks((0:fs:timescale*fs));
xticklabels({-5 : 10});
title('Z-Score Trigger 1 Heat Map')
ylabel('Trials', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition,'-','zscore colorplot'); 
title(note_title);

% save figure as matlab fig
cd (folder);
saveas(gcf,note_title, 'fig');

% mean +/- sem shadded error bar graph - trigger 1
figure; 
shadedErrorBar([-bltimewin:1/fs:timewin],zmean1,zerror1, 'lineProps','b'); hold on;
axis([-5 inf -12 20]);
plot([0;0],[-12;20],'k','Linewidth',1);
ylabel('z score', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition,'-','zscore'); 
title(note_title);

% save figure as matlab fig
cd (folder);
saveas(gcf,note_title, 'fig');
% 
%%  Quantify changes as peak (absmax) z score value

% Data are transformed into moving average of 10ms

% Time windows = 1s bin

E=(bltimewin-2)*fs;
F=(bltimewin-1)*fs; 
G=bltimewin*fs; 
H=(bltimewin+1)*fs;
I=(bltimewin+2)*fs;
J=(bltimewin+3)*fs;
K=(bltimewin+4)*fs;
L=(bltimewin+5)*fs;
M=(bltimewin+6)*fs;
N=(bltimewin+7)*fs;

%TRIGGER 1  - animal 1
zall1_1 = zeros(size(dF_trials1_1));
for i = 1:size(dF_trials1_1,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials1_1(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials1_1(i,ind)); % baseline period stdev
    zall1_1(i,:)=(dF_trials1_1(i,:) - zb)/zsd; % Z score per bin
end
zerror1_1 = std(zall1_1)/sqrt(size(zall1_1,1));
zmean1_1a =  mean (zall1_1);

%option 1: moving average option
zmean1_1 = movmean(zmean1_1a, 10);

% %option 2: binning per 10 +/- 1 ms option
% start= 1; %start at value 1
% for i= 1:((length(zmean1_1a))/10); %loops until you have 1/10th of values
%   zmean1_1(i)= mean(zmean1_1a(1, start:start+11)); %take average 
%   if start== 1
%       start= start+8;
%   else 
%       start= start+10;
%   end
% end
% 
zmean1_1A= zmean1_1 (:,E:F);
zmean1_1B= zmean1_1 (:,F:G);
zmean1_1C= zmean1_1 (:,G:H);
zmean1_1D= zmean1_1 (:,H:I);
zmean1_1E= zmean1_1 (:,I:J);
zmean1_1F= zmean1_1 (:,J:K);
zmean1_1G= zmean1_1 (:,K:L);
zmean1_1H= zmean1_1 (:,L:M);
zmean1_1I= zmean1_1 (:,M:N);

[~,X] = max(abs(zmean1_1A));
zmean1_1Apeak = zmean1_1A(X);
[~,X] = max(abs(zmean1_1B));
zmean1_1Bpeak = zmean1_1B(X);
[~,X] = max(abs(zmean1_1C));
zmean1_1Cpeak = zmean1_1C(X);
[~,X] = max(abs(zmean1_1D));
zmean1_1Dpeak = zmean1_1D(X);
[~,X] = max(abs(zmean1_1E));
zmean1_1Epeak = zmean1_1E(X);
[~,X] = max(abs(zmean1_1F));
zmean1_1Fpeak = zmean1_1F(X);
[~,X] = max(abs(zmean1_1G));
zmean1_1Gpeak = zmean1_1G(X);
[~,X] = max(abs(zmean1_1H));
zmean1_1Hpeak = zmean1_1H(X);
[~,X] = max(abs(zmean1_1I));
zmean1_1Ipeak = zmean1_1I(X);

%TRIGGER 1  - animal 2
zall1_2 = zeros(size(dF_trials1_2));
for i = 1:size(dF_trials1_2,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials1_2(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials1_2(i,ind)); % baseline period stdev
    zall1_2(i,:)=(dF_trials1_2(i,:) - zb)/zsd; % Z score per bin
end
zerror1_2 = std(zall1_2)/sqrt(size(zall1_2,1));
zmean1_2a =  mean (zall1_2);
zmean1_2 = movmean(zmean1_2a, 10);

zmean1_2A= zmean1_2 (:,E:F);
zmean1_2B= zmean1_2 (:,F:G);
zmean1_2C= zmean1_2 (:,G:H);
zmean1_2D= zmean1_2 (:,H:I);
zmean1_2E= zmean1_2 (:,I:J);
zmean1_2F= zmean1_2 (:,J:K);
zmean1_2G= zmean1_2 (:,K:L);
zmean1_2H= zmean1_2 (:,L:M);
zmean1_2I= zmean1_2 (:,M:N);

[~,X] = max(abs(zmean1_2A));
zmean1_2Apeak = zmean1_2A(X);
[~,X] = max(abs(zmean1_2B));
zmean1_2Bpeak = zmean1_2B(X);
[~,X] = max(abs(zmean1_2C));
zmean1_2Cpeak = zmean1_2C(X);
[~,X] = max(abs(zmean1_2D));
zmean1_2Dpeak = zmean1_2D(X);
[~,X] = max(abs(zmean1_2E));
zmean1_2Epeak = zmean1_2E(X);
[~,X] = max(abs(zmean1_2F));
zmean1_2Fpeak = zmean1_2F(X);
[~,X] = max(abs(zmean1_2G));
zmean1_2Gpeak = zmean1_2G(X);
[~,X] = max(abs(zmean1_2H));
zmean1_2Hpeak = zmean1_2H(X);
[~,X] = max(abs(zmean1_2I));
zmean1_2Ipeak = zmean1_2I(X);


%TRIGGER 1  - animal 3
zall1_3 = zeros(size(dF_trials1_3));
for i = 1:size(dF_trials1_3,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials1_3(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials1_3(i,ind)); % baseline period stdev
    zall1_3(i,:)=(dF_trials1_3(i,:) - zb)/zsd; % Z score per bin
end
zerror1_3 = std(zall1_3)/sqrt(size(zall1_3,1));
zmean1_3a =  mean (zall1_3);
zmean1_3 = movmean(zmean1_3a, 10);

zmean1_3A= zmean1_3 (:,E:F);
zmean1_3B= zmean1_3 (:,F:G);
zmean1_3C= zmean1_3 (:,G:H);
zmean1_3D= zmean1_3 (:,H:I);
zmean1_3E= zmean1_3 (:,I:J);
zmean1_3F= zmean1_3 (:,J:K);
zmean1_3G= zmean1_3 (:,K:L);
zmean1_3H= zmean1_3 (:,L:M);
zmean1_3I= zmean1_3 (:,M:N);

[~,X] = max(abs(zmean1_3A));
zmean1_3Apeak = zmean1_3A(X);
[~,X] = max(abs(zmean1_3B));
zmean1_3Bpeak = zmean1_3B(X);
[~,X] = max(abs(zmean1_3C));
zmean1_3Cpeak = zmean1_3C(X);
[~,X] = max(abs(zmean1_3D));
zmean1_3Dpeak = zmean1_3D(X);
[~,X] = max(abs(zmean1_3E));
zmean1_3Epeak = zmean1_3E(X);
[~,X] = max(abs(zmean1_3F));
zmean1_3Fpeak = zmean1_3F(X);
[~,X] = max(abs(zmean1_3G));
zmean1_3Gpeak = zmean1_3G(X);
[~,X] = max(abs(zmean1_3H));
zmean1_3Hpeak = zmean1_3H(X);
[~,X] = max(abs(zmean1_3I));
zmean1_3Ipeak = zmean1_3I(X);


%TRIGGER 1  - animal 4
zall1_4 = zeros(size(dF_trials1_4));
for i = 1:size(dF_trials1_4,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials1_4(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials1_4(i,ind)); % baseline period stdev
    zall1_4(i,:)=(dF_trials1_4(i,:) - zb)/zsd; % Z score per bin
end
zerror1_4 = std(zall1_4)/sqrt(size(zall1_4,1));
zmean1_4a =  mean (zall1_4);
zmean1_4 = movmean(zmean1_4a, 10);

zmean1_4A= zmean1_4 (:,E:F);
zmean1_4B= zmean1_4 (:,F:G);
zmean1_4C= zmean1_4 (:,G:H);
zmean1_4D= zmean1_4 (:,H:I);
zmean1_4E= zmean1_4 (:,I:J);
zmean1_4F= zmean1_4 (:,J:K);
zmean1_4G= zmean1_4 (:,K:L);
zmean1_4H= zmean1_4 (:,L:M);
zmean1_4I= zmean1_4 (:,M:N);

[~,X] = max(abs(zmean1_4A));
zmean1_4Apeak = zmean1_4A(X);
[~,X] = max(abs(zmean1_4B));
zmean1_4Bpeak = zmean1_4B(X);
[~,X] = max(abs(zmean1_4C));
zmean1_4Cpeak = zmean1_4C(X);
[~,X] = max(abs(zmean1_4D));
zmean1_4Dpeak = zmean1_4D(X);
[~,X] = max(abs(zmean1_4E));
zmean1_4Epeak = zmean1_4E(X);
[~,X] = max(abs(zmean1_4F));
zmean1_4Fpeak = zmean1_4F(X);
[~,X] = max(abs(zmean1_4G));
zmean1_4Gpeak = zmean1_4G(X);
[~,X] = max(abs(zmean1_4H));
zmean1_4Hpeak = zmean1_4H(X);
[~,X] = max(abs(zmean1_4I));
zmean1_4Ipeak = zmean1_4I(X);


% %TRIGGER 1  - animal 5
zall1_5 = zeros(size(dF_trials1_5));
for i = 1:size(dF_trials1_5,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials1_5(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials1_5(i,ind)); % baseline period stdev
    zall1_5(i,:)=(dF_trials1_5(i,:) - zb)/zsd; % Z score per bin
end
zerror1_5 = std(zall1_5)/sqrt(size(zall1_5,1));
zmean1_5a =  mean (zall1_5);
zmean1_5 = movmean(zmean1_5a, 10);
% 
zmean1_5A= zmean1_5 (:,E:F);
zmean1_5B= zmean1_5 (:,F:G);
zmean1_5C= zmean1_5 (:,G:H);
zmean1_5D= zmean1_5 (:,H:I);
zmean1_5E= zmean1_5 (:,I:J);
zmean1_5F= zmean1_5 (:,J:K);
zmean1_5G= zmean1_5 (:,K:L);
zmean1_5H= zmean1_5 (:,L:M);
zmean1_5I= zmean1_5 (:,M:N);

[~,X] = max(abs(zmean1_5A));
zmean1_5Apeak = zmean1_5A(X);
[~,X] = max(abs(zmean1_5B));
zmean1_5Bpeak = zmean1_5B(X);
[~,X] = max(abs(zmean1_5C));
zmean1_5Cpeak = zmean1_5C(X);
[~,X] = max(abs(zmean1_5D));
zmean1_5Dpeak = zmean1_5D(X);
[~,X] = max(abs(zmean1_5E));
zmean1_5Epeak = zmean1_5E(X);
[~,X] = max(abs(zmean1_5F));
zmean1_5Fpeak = zmean1_5F(X);
[~,X] = max(abs(zmean1_5G));
zmean1_5Gpeak = zmean1_5G(X);
[~,X] = max(abs(zmean1_5H));
zmean1_5Hpeak = zmean1_5H(X);
[~,X] = max(abs(zmean1_5I));
zmean1_5Ipeak = zmean1_5I(X);

%
% %TRIGGER 1  - animal 6
zall1_6 = zeros(size(dF_trials1_6));
for i = 1:size(dF_trials1_6,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials1_6(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials1_6(i,ind)); % baseline period stdev
    zall1_6(i,:)=(dF_trials1_6(i,:) - zb)/zsd; % Z score per bin
end
zerror1_6 = std(zall1_6)/sqrt(size(zall1_6,1));
zmean1_6a =  mean (zall1_6);
zmean1_6 = movmean(zmean1_6a, 10);
% 
zmean1_6A= zmean1_6 (:,E:F);
zmean1_6B= zmean1_6 (:,F:G);
zmean1_6C= zmean1_6 (:,G:H);
zmean1_6D= zmean1_6 (:,H:I);
zmean1_6E= zmean1_6 (:,I:J);
zmean1_6F= zmean1_6 (:,J:K);
zmean1_6G= zmean1_6 (:,K:L);
zmean1_6H= zmean1_6 (:,L:M);
zmean1_6I= zmean1_6 (:,M:N);

[~,X] = max(abs(zmean1_6A));
zmean1_6Apeak = zmean1_6A(X);
[~,X] = max(abs(zmean1_6B));
zmean1_6Bpeak = zmean1_6B(X);
[~,X] = max(abs(zmean1_6C));
zmean1_6Cpeak = zmean1_6C(X);
[~,X] = max(abs(zmean1_6D));
zmean1_6Dpeak = zmean1_6D(X);
[~,X] = max(abs(zmean1_6E));
zmean1_6Epeak = zmean1_6E(X);
[~,X] = max(abs(zmean1_6F));
zmean1_6Fpeak = zmean1_6F(X);
[~,X] = max(abs(zmean1_6G));
zmean1_6Gpeak = zmean1_6G(X);
[~,X] = max(abs(zmean1_6H));
zmean1_6Hpeak = zmean1_6H(X);
[~,X] = max(abs(zmean1_6I));
zmean1_6Ipeak = zmean1_6I(X);


% TRIGGER 1  - animal 7
zall1_7 = zeros(size(dF_trials1_7));
for i = 1:size(dF_trials1_7,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials1_7(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials1_7(i,ind)); % baseline period stdev
    zall1_7(i,:)=(dF_trials1_7(i,:) - zb)/zsd; % Z score per bin
end
zerror1_7 = std(zall1_7)/sqrt(size(zall1_7,1));
zmean1_7a =  mean (zall1_7);
zmean1_7 = movmean(zmean1_7a, 10);

zmean1_7A= zmean1_7 (:,E:F);
zmean1_7B= zmean1_7 (:,F:G);
zmean1_7C= zmean1_7 (:,G:H);
zmean1_7D= zmean1_7 (:,H:I);
zmean1_7E= zmean1_7 (:,I:J);
zmean1_7F= zmean1_7 (:,J:K);
zmean1_7G= zmean1_7 (:,K:L);
zmean1_7H= zmean1_7 (:,L:M);
zmean1_7I= zmean1_7 (:,M:N);

[~,X] = max(abs(zmean1_7A));
zmean1_7Apeak = zmean1_7A(X);
[~,X] = max(abs(zmean1_7B));
zmean1_7Bpeak = zmean1_7B(X);
[~,X] = max(abs(zmean1_7C));
zmean1_7Cpeak = zmean1_7C(X);
[~,X] = max(abs(zmean1_7D));
zmean1_7Dpeak = zmean1_7D(X);
[~,X] = max(abs(zmean1_7E));
zmean1_7Epeak = zmean1_7E(X);
[~,X] = max(abs(zmean1_7F));
zmean1_7Fpeak = zmean1_7F(X);
[~,X] = max(abs(zmean1_7G));
zmean1_7Gpeak = zmean1_7G(X);
[~,X] = max(abs(zmean1_7H));
zmean1_7Hpeak = zmean1_7H(X);
[~,X] = max(abs(zmean1_7I));
zmean1_7Ipeak = zmean1_7I(X);


%TRIGGER 1  - animal 8
zall1_8 = zeros(size(dF_trials1_8));
for i = 1:size(dF_trials1_8,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials1_8(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials1_8(i,ind)); % baseline period stdev
    zall1_8(i,:)=(dF_trials1_8(i,:) - zb)/zsd; % Z score per bin
end
zerror1_8 = std(zall1_8)/sqrt(size(zall1_8,1));
zmean1_8a =  mean (zall1_8);
zmean1_8 = movmean(zmean1_8a, 10);
% 
zmean1_8A= zmean1_8 (:,E:F);
zmean1_8B= zmean1_8 (:,F:G);
zmean1_8C= zmean1_8 (:,G:H);
zmean1_8D= zmean1_8 (:,H:I);
zmean1_8E= zmean1_8 (:,I:J);
zmean1_8F= zmean1_8 (:,J:K);
zmean1_8G= zmean1_8 (:,K:L);
zmean1_8H= zmean1_8 (:,L:M);
zmean1_8I= zmean1_8 (:,M:N);

[~,X] = max(abs(zmean1_8A));
zmean1_8Apeak = zmean1_8A(X);
[~,X] = max(abs(zmean1_8B));
zmean1_8Bpeak = zmean1_8B(X);
[~,X] = max(abs(zmean1_8C));
zmean1_8Cpeak = zmean1_8C(X);
[~,X] = max(abs(zmean1_8D));
zmean1_8Dpeak = zmean1_8D(X);
[~,X] = max(abs(zmean1_8E));
zmean1_8Epeak = zmean1_8E(X);
[~,X] = max(abs(zmean1_8F));
zmean1_8Fpeak = zmean1_8F(X);
[~,X] = max(abs(zmean1_8G));
zmean1_8Gpeak = zmean1_8G(X);
[~,X] = max(abs(zmean1_8H));
zmean1_8Hpeak = zmean1_8H(X);
[~,X] = max(abs(zmean1_8I));
zmean1_8Ipeak = zmean1_8I(X);


% 
zabsmax_PreStim21= [zmean1_1Apeak;zmean1_2Apeak;zmean1_3Apeak;zmean1_4Apeak;zmean1_5Apeak;zmean1_6Apeak;zmean1_7Apeak;zmean1_8Apeak];
zabsmax_Prestim10= [zmean1_1Bpeak;zmean1_2Bpeak;zmean1_3Bpeak;zmean1_4Bpeak;zmean1_5Bpeak;zmean1_6Bpeak;zmean1_7Bpeak;zmean1_8Bpeak];
zabsmax_PostStim01= [zmean1_1Cpeak;zmean1_2Cpeak;zmean1_3Cpeak;zmean1_4Cpeak;zmean1_5Cpeak;zmean1_6Cpeak;zmean1_7Cpeak;zmean1_8Cpeak];
zabsmax_PostStim12= [zmean1_1Dpeak;zmean1_2Dpeak;zmean1_3Dpeak;zmean1_4Dpeak;zmean1_5Dpeak;zmean1_6Dpeak;zmean1_7Dpeak;zmean1_8Dpeak];
zabsmax_PostStim23= [zmean1_1Epeak;zmean1_2Epeak;zmean1_3Epeak;zmean1_4Epeak;zmean1_5Epeak;zmean1_6Epeak;zmean1_7Epeak;zmean1_8Epeak];
zabsmax_PostStim34= [zmean1_1Fpeak;zmean1_2Fpeak;zmean1_3Fpeak;zmean1_4Fpeak;zmean1_5Fpeak;zmean1_6Fpeak;zmean1_7Fpeak;zmean1_8Fpeak];
zabsmax_PostStim45= [zmean1_1Gpeak;zmean1_2Gpeak;zmean1_3Gpeak;zmean1_4Gpeak;zmean1_5Gpeak;zmean1_6Gpeak;zmean1_7Gpeak;zmean1_8Gpeak];
zabsmax_PostStim56= [zmean1_1Hpeak;zmean1_2Hpeak;zmean1_3Hpeak;zmean1_4Hpeak;zmean1_5Hpeak;zmean1_6Hpeak;zmean1_7Hpeak;zmean1_8Hpeak];
zabsmax_PostStim67= [zmean1_1Ipeak;zmean1_2Ipeak;zmean1_3Ipeak;zmean1_4Ipeak;zmean1_5Ipeak;zmean1_6Ipeak;zmean1_7Ipeak;zmean1_8Ipeak];

 

%% Save matlab data file
% You can save space by saving only necessary data - deleting unwanted data in the
% workspace - before saving .mat file

% cd 'filepath';
% save filename.mat
