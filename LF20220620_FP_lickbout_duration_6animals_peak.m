%% Script written by Lauren Faget - Hnasko lab - July 2021
% For Fiber photometry data analysis

% PLEASE DO NOT SHARE OUTSIDE OF THE HNASKO LAB

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
note_condition1 = 'bouts of licks - onset';
note_condition2 = 'bouts of licks - offset';

%% MODIFY trigger names here
%trigger names can be found in data/ epocs/ 
%format is '3 letters - underscore' or 'x - 3 numbers/numbers-letters combo - underscore'

EPOC_1 = 'RRli'; % enter name of trigger 1
EPOC_2 = 'LLli'; % enter name of trigger 1
% EPOC_3 = 'Rbo_'; % enter name of trigger 1
% EPOC_4 = 'LLbo'; % enter name of trigger 1

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
trt1_1=(tr1_1*fs);
trt1_2=(tr1_2*fs);
trt1_3=(tr1_3*fs);
trt1_4=(tr1_4*fs);
trt1_5=(tr1_5*fs);
trt1_6=(tr1_6*fs);

%
%%  Extract bout onset/offset/duration from lick epoc arrays

%animal 1
i=1;
for i= 1:length (tr1_1)-1
    if (tr1_1(i+1)-tr1_1(i)) <= 0.75
    newmatrix (i)= tr1_1(i);
    else newmatrix (i)= 0;
    end
 i= i+1;    
end

i=1;
for i= 1:length (newmatrix)-1
    if (newmatrix(i+1) - newmatrix(i)) == newmatrix(i+1) && (newmatrix(i+1) - newmatrix(i)) ~=0
    data1.epocs.realbout.onset (i+1) = newmatrix (i+1);
    else data1.epocs.realbout.onset (i+1)= NaN;
    end
 i= i+1;
end
   if newmatrix(1)==0
   data1.epocs.realbout.onset(1)= NaN;
   else data1.epocs.realbout.onset(1)= newmatrix(1);
   end
   
data1.epocs.realbout.onset(data1.epocs.realbout.onset==0) = NaN;
data1.epocs.realbout.onset = rmmissing(data1.epocs.realbout.onset);
data1.epocs.realbout.onset = data1.epocs.realbout.onset (1:end-1);

i=1;
for i= 1:length (newmatrix)-1
    if (newmatrix(i+1) - newmatrix(i)) == -newmatrix(i) && newmatrix(i) ~= 0
    data1.epocs.realbout.offset (i) = tr1_1 (i+1);
    else data1.epocs.realbout.offset (i)= NaN;
    end
 i= i+1;
end
data1.epocs.realbout.offset(data1.epocs.realbout.offset==0) = NaN;
data1.epocs.realbout.offset = rmmissing(data1.epocs.realbout.offset);

%animal 2
i=1;
for i= 1:length (tr1_2)-1
    if (tr1_2(i+1)-tr1_2(i)) <= 0.75
    newmatrix2 (i)= tr1_2(i);
    else newmatrix2 (i)= 0;
    end
 i= i+1;    
end

i=1;
for i= 1:length (newmatrix2)-1
    if (newmatrix2(i+1) - newmatrix2(i)) == newmatrix2(i+1) && (newmatrix2(i+1) - newmatrix2(i)) ~=0
    data2.epocs.realbout.onset (i+1) = newmatrix2 (i+1);
    else data2.epocs.realbout.onset (i+1)= NaN;
    end
 i= i+1;
end
   if newmatrix2(1)==0
   data2.epocs.realbout.onset(1)= NaN;
   else data2.epocs.realbout.onset(1)= newmatrix2(1);
   end

data2.epocs.realbout.onset(data2.epocs.realbout.onset==0) = NaN;
data2.epocs.realbout.onset = rmmissing(data2.epocs.realbout.onset);
data2.epocs.realbout.onset = data2.epocs.realbout.onset (1:end-1);

i=1;
for i= 1:length (newmatrix2)-1
    if (newmatrix2(i+1) - newmatrix2(i)) == -newmatrix2(i) && newmatrix2(i) ~= 0
    data2.epocs.realbout.offset (i) = tr1_2 (i+1);
    else data2.epocs.realbout.offset (i)= NaN;
    end
 i= i+1;
end
data2.epocs.realbout.offset(data2.epocs.realbout.offset==0) = NaN;
data2.epocs.realbout.offset = rmmissing(data2.epocs.realbout.offset);

%animal 3
i=1;
for i= 1:length (tr1_3)-1
    if (tr1_3(i+1)-tr1_3(i)) <= 0.75
    newmatrix3 (i)= tr1_3(i);
    else newmatrix3 (i)= 0;
    end
 i= i+1;    
end

i=1;
for i= 1:length (newmatrix3)-1
    if (newmatrix3(i+1) - newmatrix3(i)) == newmatrix3(i+1) && (newmatrix3(i+1) - newmatrix3(i)) ~=0
    data3.epocs.realbout.onset (i+1) = newmatrix3 (i+1);
    else data3.epocs.realbout.onset (i+1)= NaN;
    end
 i= i+1;
end
   if newmatrix3(1)==0
   data3.epocs.realbout.onset(1)= NaN;
   else data3.epocs.realbout.onset(1)= newmatrix3(1);
   end

data3.epocs.realbout.onset(data3.epocs.realbout.onset==0) = NaN;
data3.epocs.realbout.onset = rmmissing(data3.epocs.realbout.onset);
data3.epocs.realbout.onset = data3.epocs.realbout.onset (1:end-1);

i=1;
for i= 1:length (newmatrix3)-1
    if (newmatrix3(i+1) - newmatrix3(i)) == -newmatrix3(i) && newmatrix3(i) ~= 0
    data3.epocs.realbout.offset (i) = tr1_3 (i+1);
    else data3.epocs.realbout.offset (i)= NaN;
    end
 i= i+1;
end
data3.epocs.realbout.offset(data3.epocs.realbout.offset==0) = NaN;
data3.epocs.realbout.offset = rmmissing(data3.epocs.realbout.offset);

%animal4
i=1;
for i= 1:length (tr1_4)-1
    if (tr1_4(i+1)-tr1_4(i)) <= 0.75
    newmatrix4 (i)= tr1_4(i);
    else newmatrix4 (i)= 0;
    end
 i= i+1;    
end

i=1;
for i= 1:length (newmatrix4)-1
    if (newmatrix4(i+1) - newmatrix4(i)) == newmatrix4(i+1) && (newmatrix4(i+1) - newmatrix4(i)) ~=0
    data4.epocs.realbout.onset (i+1) = newmatrix4 (i+1);
    else data4.epocs.realbout.onset (i+1)= NaN;
    end
 i= i+1;
end
   if newmatrix4(1)==0
   data4.epocs.realbout.onset(1)= NaN;
   else data4.epocs.realbout.onset(1)= newmatrix4(1);
   end

data4.epocs.realbout.onset(data4.epocs.realbout.onset==0) = NaN;
data4.epocs.realbout.onset = rmmissing(data4.epocs.realbout.onset);
data4.epocs.realbout.onset = data4.epocs.realbout.onset (1:end-1);

i=1;
for i= 1:length (newmatrix4)-1
    if (newmatrix4(i+1) - newmatrix4(i)) == -newmatrix4(i) && newmatrix4(i) ~= 0
    data4.epocs.realbout.offset (i) = tr1_4 (i+1);
    else data4.epocs.realbout.offset (i)= NaN;
    end
 i= i+1;
end
data4.epocs.realbout.offset(data4.epocs.realbout.offset==0) = NaN;
data4.epocs.realbout.offset = rmmissing(data4.epocs.realbout.offset);

%animal 5
i=1;
for i= 1:length (tr1_5)-1
    if (tr1_5(i+1)-tr1_5(i)) <= 0.75
    newmatrix5 (i)= tr1_5(i);
    else newmatrix5 (i)= 0;
    end
 i= i+1;    
end

i=1;
for i= 1:length (newmatrix5)-1
    if (newmatrix5(i+1) - newmatrix5(i)) == newmatrix5(i+1) && (newmatrix5(i+1) - newmatrix5(i)) ~=0
    data5.epocs.realbout.onset (i+1) = newmatrix5 (i+1);
    else data5.epocs.realbout.onset (i+1)= NaN;
    end
 i= i+1;
end
   if newmatrix5(1)==0
   data5.epocs.realbout.onset(1)= NaN;
   else data5.epocs.realbout.onset(1)= newmatrix5(1);
   end

data5.epocs.realbout.onset(data5.epocs.realbout.onset==0) = NaN;
data5.epocs.realbout.onset = rmmissing(data5.epocs.realbout.onset);
data5.epocs.realbout.onset = data5.epocs.realbout.onset (1:end-1);

i=1;
for i= 1:length (newmatrix5)-1
    if (newmatrix5(i+1) - newmatrix5(i)) == -newmatrix5(i) && newmatrix5(i) ~= 0
    data5.epocs.realbout.offset (i) = tr1_5 (i+1);
    else data5.epocs.realbout.offset (i)= NaN;
    end
 i= i+1;
end
data5.epocs.realbout.offset(data5.epocs.realbout.offset==0) = NaN;
data5.epocs.realbout.offset = rmmissing(data5.epocs.realbout.offset);

%animal 6
i=1;
for i= 1:length (tr1_6)-1
    if (tr1_6(i+1)-tr1_6(i)) <= 0.75
    newmatrix6 (i)= tr1_6(i);
    else newmatrix6 (i)= 0;
    end
 i= i+1;    
end

i=1;
for i= 1:length (newmatrix6)-1
    if (newmatrix6(i+1) - newmatrix6(i)) == newmatrix6(i+1) && (newmatrix6(i+1) - newmatrix6(i)) ~=0
    data6.epocs.realbout.onset (i+1) = newmatrix6 (i+1);
    else data6.epocs.realbout.onset (i+1)= NaN;
    end
 i= i+1;
end
   if newmatrix6(1)==0
   data6.epocs.realbout.onset(1)= NaN;
   else data6.epocs.realbout.onset(1)= newmatrix6(1);
   end

data6.epocs.realbout.onset(data6.epocs.realbout.onset==0) = NaN;
data6.epocs.realbout.onset = rmmissing(data6.epocs.realbout.onset);
data6.epocs.realbout.onset = data6.epocs.realbout.onset (1:end-1);

i=1;
for i= 1:length (newmatrix6)-1
    if (newmatrix6(i+1) - newmatrix6(i)) == -newmatrix6(i) && newmatrix6(i) ~= 0
    data6.epocs.realbout.offset (i) = tr1_6 (i+1);
    else data6.epocs.realbout.offset (i)= NaN;
    end
 i= i+1;
end
data6.epocs.realbout.offset(data6.epocs.realbout.offset==0) = NaN;
data6.epocs.realbout.offset = rmmissing(data6.epocs.realbout.offset);

%
tr3_1=data1.epocs.realbout.onset';
tr3_2=data2.epocs.realbout.onset';
tr3_3=data3.epocs.realbout.onset';
tr3_4=data4.epocs.realbout.onset';
tr3_5=data5.epocs.realbout.onset';
tr3_6=data6.epocs.realbout.onset';
% 
trt3_1=(tr3_1*fs);
trt3_2=(tr3_2*fs);
trt3_3=(tr3_3*fs);
trt3_4=(tr3_4*fs);
trt3_5=(tr3_5*fs);
trt3_6=(tr3_6*fs);

tr4_1=data1.epocs.realbout.offset';
tr4_2=data2.epocs.realbout.offset';
tr4_3=data3.epocs.realbout.offset';
tr4_4=data4.epocs.realbout.offset';
tr4_5=data5.epocs.realbout.offset';
tr4_6=data6.epocs.realbout.offset';
% 
trt4_1=(tr4_1*fs);
trt4_2=(tr4_2*fs);
trt4_3=(tr4_3*fs);
trt4_4=(tr4_4*fs);
trt4_5=(tr4_5*fs);
trt4_6=(tr4_6*fs);


% animal 1
i=1;
for i= 1:length (tr3_1)
    data1.epocs.realbout.data (i)= data1.epocs.realbout.offset (i)- data1.epocs.realbout.onset (i);
 i= i+1;
end
%
% animal 2
i=1;
for i= 1:length (tr3_2)
    data2.epocs.realbout.data (i)= data2.epocs.realbout.offset (i)- data2.epocs.realbout.onset (i);
 i= i+1;
end
%
% animal 3
i=1;
for i= 1:length (tr3_3)
    data3.epocs.realbout.data (i)= data3.epocs.realbout.offset (i)- data3.epocs.realbout.onset (i);
 i= i+1;
end
%
% animal 4
i=1;
for i= 1:length (tr3_4)
    data4.epocs.realbout.data (i)= data4.epocs.realbout.offset (i)- data4.epocs.realbout.onset (i);
 i= i+1;
end
%
% animal 5
i=1;
for i= 1:length (tr3_5)
    data5.epocs.realbout.data (i)= data5.epocs.realbout.offset (i)- data5.epocs.realbout.onset (i);
 i= i+1;
end
%
% animal 6
i=1;
for i= 1:length (tr3_6)
    data6.epocs.realbout.data (i)= data6.epocs.realbout.offset (i)- data6.epocs.realbout.onset (i);
 i= i+1;
end


tr3_1length=data1.epocs.realbout.data';
tr3_2length=data2.epocs.realbout.data';
tr3_3length=data3.epocs.realbout.data';
tr3_4length=data4.epocs.realbout.data';
tr3_5length=data5.epocs.realbout.data';
tr3_6length=data6.epocs.realbout.data';

%% Plot all licks, bout onset & bout offset

figure;
title('Entire Experiment - animal 1');hold on;
plot((1/fs:1/fs:length(dF1)/fs),dF1,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_1)
     plot([tr1_1(i);tr1_1(i)],[-15;15],'c');
 end
  for i = 1:length(tr3_1)
     plot([tr3_1(i);tr3_1(i)],[-15;15],'r');
 end
  for i = 1:length(tr4_1)
     plot([tr4_1(i);tr4_1(i)],[-15;15],'k');
 end

 figure;
title('Entire Experiment - animal 2');hold on;
plot((1/fs:1/fs:length(dF2)/fs),dF2,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

 for i = 1:length(tr1_2)
     plot([tr1_2(i);tr1_2(i)],[-15;15],'c');
 end
  for i = 1:length(tr3_2)
     plot([tr3_2(i);tr3_2(i)],[-15;15],'r');
 end
  for i = 1:length(tr4_2)
     plot([tr4_2(i);tr4_2(i)],[-15;15],'k');
 end

figure;
title('Entire Experiment - animal 3');hold on;
plot((1/fs:1/fs:length(dF3)/fs),dF3,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_3)
     plot([tr1_3(i);tr1_3(i)],[-15;15],'c');
 end
  for i = 1:length(tr3_3)
     plot([tr3_3(i);tr3_3(i)],[-15;15],'r');
 end
  for i = 1:length(tr4_3)
     plot([tr4_3(i);tr4_3(i)],[-15;15],'k');
 end

figure;
title('Entire Experiment - animal 4');hold on;
plot((1/fs:1/fs:length(dF4)/fs),dF4,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i = 1:length(tr1_4)
     plot([tr1_4(i);tr1_4(i)],[-15;15],'c');
 end
  for i = 1:length(tr3_4)
     plot([tr3_4(i);tr3_4(i)],[-15;15],'r');
 end
  for i = 1:length(tr4_4)
     plot([tr4_4(i);tr4_4(i)],[-15;15],'k');
 end

figure;
title('Entire Experiment - animal 5');hold on;
plot((1/fs:1/fs:length(dF5)/fs),dF5,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;
  for i = 1:length(tr1_5)
     plot([tr1_5(i);tr1_5(i)],[-15;15],'c');
 end
  for i = 1:length(tr3_5)
     plot([tr3_5(i);tr3_5(i)],[-15;15],'r');
 end
  for i = 1:length(tr4_5)
     plot([tr4_5(i);tr4_5(i)],[-15;15],'k');
 end

 figure;
title('Entire Experiment - animal 6');hold on;
plot((1/fs:1/fs:length(dF6)/fs),dF6,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;
  for i = 1:length(tr1_6)
     plot([tr1_6(i);tr1_6(i)],[-15;15],'c');
 end
  for i = 1:length(tr3_6)
     plot([tr3_6(i);tr3_6(i)],[-15;15],'r');
 end
  for i = 1:length(tr4_6)
     plot([tr4_6(i);tr4_6(i)],[-15;15],'k');
 end

% close all;

%% C.this plots +/- time windows around bout onset & offset

%DEFINE your time window here;
bltimewin=5;  %% baseline time window in seconds
timewin=10;   %% time window after trigger in seconds
timescale = bltimewin + timewin;

dF_trials3_1=[]; 
dF_trials3_2=[]; 
dF_trials3_3=[]; 
dF_trials3_4=[]; 
dF_trials3_5=[]; 
dF_trials3_6=[]; 

dF_trials4_1=[]; 
dF_trials4_2=[]; 
dF_trials4_3=[]; 
dF_trials4_4=[]; 
dF_trials4_5=[]; 
dF_trials4_6=[]; 

%animal 1
% extract trials for the third trigger.
count=0;
for i = 1:length(trt3_1)
     if (trt3_1(i)+timewin*fs<length(dF1) && trt3_1(i)-timewin*fs>0)     
        count=count+1;
        dF_trials3_1(count,:)=dF1(trt3_1(i)-(bltimewin*fs):trt3_1(i)+(timewin*fs));  
    end
end
% extract trials for the fourth trigger.
count=0;
for i = 1:length(trt4_1)
     if (trt4_1(i)+timewin*fs<length(dF1) && trt4_1(i)-timewin*fs>0)     
        count=count+1;
        dF_trials4_1(count,:)=dF1(trt4_1(i)-(bltimewin*fs):trt4_1(i)+(timewin*fs));  
    end
end

%animal 2
% extract trials for the third trigger.
count=0;
for i = 1:length(trt3_2)
     if (trt3_2(i)+timewin*fs<length(dF2) && trt3_2(i)-timewin*fs>0)     
        count=count+1;
        dF_trials3_2(count,:)=dF2(trt3_2(i)-(bltimewin*fs):trt3_2(i)+(timewin*fs));  
    end
end
% extract trials for the fourth trigger.
count=0;
for i = 1:length(trt4_2)
     if (trt4_2(i)+timewin*fs<length(dF2) && trt4_2(i)-timewin*fs>0)     
        count=count+1;
        dF_trials4_2(count,:)=dF2(trt4_2(i)-(bltimewin*fs):trt4_2(i)+(timewin*fs));  
    end
end

%animal 3
% extract trials for the third trigger.
count=0;
for i = 1:length(trt3_3)
     if (trt3_3(i)+timewin*fs<length(dF3) && trt3_3(i)-timewin*fs>0)     
        count=count+1;
        dF_trials3_3(count,:)=dF3(trt3_3(i)-(bltimewin*fs):trt3_3(i)+(timewin*fs));  
    end
end
% extract trials for the fourth trigger.
count=0;
for i = 1:length(trt4_3)
     if (trt4_3(i)+timewin*fs<length(dF3) && trt4_3(i)-timewin*fs>0)     
        count=count+1;
        dF_trials4_3(count,:)=dF3(trt4_3(i)-(bltimewin*fs):trt4_3(i)+(timewin*fs));  
    end
end

%animal 4
% extract trials for the third trigger.
count=0;
for i = 1:length(trt3_4)
     if (trt3_4(i)+timewin*fs<length(dF4) && trt3_4(i)-timewin*fs>0)     
        count=count+1;
        dF_trials3_4(count,:)=dF4(trt3_4(i)-(bltimewin*fs):trt3_4(i)+(timewin*fs));  
    end
end
% extract trials for the fourth trigger.
count=0;
for i = 1:length(trt4_4)
     if (trt4_4(i)+timewin*fs<length(dF4) && trt4_4(i)-timewin*fs>0)     
        count=count+1;
        dF_trials4_4(count,:)=dF4(trt4_4(i)-(bltimewin*fs):trt4_4(i)+(timewin*fs));  
    end
end

%animal 5
% extract trials for the third trigger.
count=0;
for i = 1:length(trt3_5)
     if (trt3_5(i)+timewin*fs<length(dF5) && trt3_5(i)-timewin*fs>0)     
        count=count+1;
        dF_trials3_5(count,:)=dF5(trt3_5(i)-(bltimewin*fs):trt3_5(i)+(timewin*fs));  
    end
end
% extract trials for the fourth trigger.
count=0;
for i = 1:length(trt4_5)
     if (trt4_5(i)+timewin*fs<length(dF5) && trt4_5(i)-timewin*fs>0)     
        count=count+1;
        dF_trials4_5(count,:)=dF5(trt4_5(i)-(bltimewin*fs):trt4_5(i)+(timewin*fs));  
    end
end

%animal 6
% extract trials for the third trigger.
count=0;
for i = 1:length(trt3_6)
     if (trt3_6(i)+timewin*fs<length(dF6) && trt3_6(i)-timewin*fs>0)     
        count=count+1;
        dF_trials3_6(count,:)=dF6(trt3_6(i)-(bltimewin*fs):trt3_6(i)+(timewin*fs));  
    end
end
% extract trials for the fourth trigger.
count=0;
for i = 1:length(trt4_6)
     if (trt4_6(i)+timewin*fs<length(dF6) && trt4_6(i)-timewin*fs>0)     
        count=count+1;
        dF_trials4_6(count,:)=dF6(trt4_6(i)-(bltimewin*fs):trt4_6(i)+(timewin*fs));  
    end
end


%Non-normalized dF/F for all trials and all animals
dF_trials3= [dF_trials3_1;dF_trials3_2;dF_trials3_3;dF_trials3_4;dF_trials3_5;dF_trials3_6];
dF_trials4= [dF_trials4_1;dF_trials4_2;dF_trials4_3;dF_trials4_4;dF_trials4_5;dF_trials4_6];


%% 
% %%  extract Normalized data per trigger and graph ColorBars
% 
% %VALUES TO MODIFY ACCORDING TO NEED HERE
% A=(bltimewin-1)*fs; % beginning of baseline to normalize with
% B=bltimewin*fs; % end of baseline to normalize with
% % time unit is in ms. Multiply value by fs to obtain seconds.
% % normalize data using a baseline 1 second prior to trigger in this
% % example
% 
% BASELINE_NORM = [A B];
% 
% bl3_1=dF_trials3_1(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
% blm3_1=mean(bl3_1,2);
% dFn3_1=(dF_trials3_1-repmat(blm3_1,1,size(dF_trials3_1,2)));
% 
% bl3_2=dF_trials3_2(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
% blm3_2=mean(bl3_2,2);
% dFn3_2=(dF_trials3_2-repmat(blm3_2,1,size(dF_trials3_2,2)));
% 
% bl3_3=dF_trials3_3(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
% blm3_3=mean(bl3_3,2);
% dFn3_3=(dF_trials3_3-repmat(blm3_3,1,size(dF_trials3_3,2)));
% 
% bl3_4=dF_trials3_4(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
% blm3_4=mean(bl3_4,2);
% dFn3_4=(dF_trials3_4-repmat(blm3_4,1,size(dF_trials3_4,2)));
% 
% bl3_5=dF_trials3_5(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
% blm3_5=mean(bl3_5,2);
% dFn3_5=(dF_trials3_5-repmat(blm3_5,1,size(dF_trials3_5,2)));
% 
% bl3_6=dF_trials3_6(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
% blm3_6=mean(bl3_6,2);
% dFn3_6=(dF_trials3_6-repmat(blm3_6,1,size(dF_trials3_6,2)));
% 
% 
% %Normalized dF/F for all trials and all animals
% dFn3= [dFn3_1;dFn3_2;dFn3_3;dFn3_4;dFn3_5;dFn3_6];
% 
% 
% figure;
% imagesc(dFn3); colorbar; hold on;
% caxis([-4 8]);  % Y axis scale
% xline(5*fs,'-k','linewidth',2); % create a line at trigger onset
% xticks((0:fs:timescale*fs));  % X axis scale
% xticklabels({-5 : 10});  % X axis scale labeling
% title('normalized dF/F trigger 1');
% ylabel('Trials', 'FontSize', 12);
% xlabel('Time (s)','FontSize',12)
% set (gca,'TickDir','out');
% 

%% Z SCORE - bout onset

%VALUES TO MODIFY ACCORDING TO NEED HERE
C=1; % beginning of baseline  (1 means first data point of the 5s baseline)
D=bltimewin*fs; % end of baseline
% time unit is in ms. Multiply value by fs to obtain seconds.
% This example use a baseline of 5 seconds prior to trigger

BASELINE_Z = [C D]; 

%TRIGGER 1
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
note_title =  strcat(file,'-',experiment,'-',note_condition1,'-','zscore colorplot'); 
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
note_title =  strcat(file,'-',experiment,'-',note_condition1,'-','zscore'); 
title(note_title);

% save figure as matlab fig
% cd (folder);
% saveas(gcf,note_title, 'fig');
% 
%% Z SCORE - bout offset

%VALUES TO MODIFY ACCORDING TO NEED HERE
C=1; % beginning of baseline  (1 means first data point of the 5s baseline)
D=bltimewin*fs; % end of baseline
% time unit is in ms. Multiply value by fs to obtain seconds.
% This example use a baseline of 5 seconds prior to trigger

BASELINE_Z = [C D]; 

%TRIGGER 4
zall3 = zeros(size(dF_trials4));
for i = 1:size(dF_trials4,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials4(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials4(i,ind)); % baseline period stdev
    zall4(i,:)=(dF_trials4(i,:) - zb)/zsd; % Z score per bin
end
zerror4 = std(zall4)/sqrt(size(zall4,1));
zmean4 =  mean (zall4);

%FIGURES
% heat map - colorbar graph - trigger 4
figure; 
imagesc(zall4); colorbar; hold on;
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

% mean +/- sem shadded error bar graph - trigger 4
figure; 
shadedErrorBar([-bltimewin:1/fs:timewin],zmean4,zerror4, 'lineProps','b'); hold on;
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

%%  Quantify changes as min/max z score value for bout onset
% per 1s bin
C=1; % beginning of baseline  (1 means first data point of the 5s baseline)
D=bltimewin*fs; % end of baseline
% time unit is in ms. Multiply value by fs to obtain seconds.
% This example use a baseline of 5 seconds prior to trigger

BASELINE_Z = [C D]; 

E=bltimewin*fs; 
F=(bltimewin+1)*fs;

%TRIGGER 3  - animal 1
zall3_1 = zeros(size(dF_trials3_1));
for i = 1:size(dF_trials3_1,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_1(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_1(i,ind)); % baseline period stdev
    zall3_1(i,:)=(dF_trials3_1(i,:) - zb)/zsd; % Z score per bin
    zall3_1A (i,:)= zall3_1 (i,E:F); % extract 1s time window
    zall3_1Amm (i,:)= movmean(zall3_1A (i,:), 10); % transform data into moving average of 10ms
    [~,X]= max(abs(zall3_1Amm(i,:))); % find peak (positive or negative value) in that 1s time window
    zall3_1Apeak (i)= zall3_1Amm(i,X); % extract this peak value and put it in a new matrix
end
zall3_1Apeak=zall3_1Apeak';


%TRIGGER 3  - animal 2
zall3_2 = zeros(size(dF_trials3_2));
for i = 1:size(dF_trials3_2,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_2(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_2(i,ind)); % baseline period stdev
    zall3_2(i,:)=(dF_trials3_2(i,:) - zb)/zsd; % Z score per bin
    zall3_2A (i,:)= zall3_2 (i,E:F); % extract 1s time window
    zall3_2Amm (i,:)= movmean(zall3_2A(i,:), 10); % transform data into moving average of 10ms
    [~,X]= max(abs(zall3_2Amm(i,:))); % find peak (positive or negative value) in that 1s time window
    zall3_2Apeak (i)= zall3_2Amm(i,X); % extract this peak value and put it in a new matrix
end
zall3_2Apeak=zall3_2Apeak';


%TRIGGER 3  - animal 3
zall3_3 = zeros(size(dF_trials3_3));
for i = 1:size(dF_trials3_3,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_3(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_3(i,ind)); % baseline period stdev
    zall3_3(i,:)=(dF_trials3_3(i,:) - zb)/zsd; % Z score per bin
    zall3_3A (i,:)= zall3_3 (i,E:F); % extract 1s time window
    zall3_3Amm (i,:)= movmean(zall3_3A(i,:), 10); % transform data into moving average of 10ms
    [~,X]= max(abs(zall3_3Amm(i,:))); % find peak (positive or negative value) in that 1s time window
    zall3_3Apeak (i)= zall3_3Amm(i,X); % extract this peak value and put it in a new matrix
end
zall3_3Apeak=zall3_3Apeak';


%TRIGGER 3  - animal 4
zall3_4 = zeros(size(dF_trials3_4));
for i = 1:size(dF_trials3_4,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_4(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_4(i,ind)); % baseline period stdev
    zall3_4(i,:)=(dF_trials3_4(i,:) - zb)/zsd; % Z score per bin
    zall3_4A (i,:)= zall3_4 (i,E:F); % extract 1s time window
    zall3_4Amm (i,:)= movmean(zall3_4A(i,:), 10); % transform data into moving average of 10ms
    [~,X]= max(abs(zall3_4Amm(i,:))); % find peak (positive or negative value) in that 1s time window
    zall3_4Apeak (i)= zall3_4Amm(i,X); % extract this peak value and put it in a new matrix
end
zall3_4Apeak=zall3_4Apeak';


% %TRIGGER 3  - animal 5
zall3_5 = zeros(size(dF_trials3_5));
for i = 1:size(dF_trials3_5,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_5(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_5(i,ind)); % baseline period stdev
    zall3_5(i,:)=(dF_trials3_5(i,:) - zb)/zsd; % Z score per bin
    zall3_5A (i,:)= zall3_5 (i,E:F); % extract 1s time window
    zall3_5Amm (i,:)= movmean(zall3_5A(i,:), 10); % transform data into moving average of 10ms
    [~,X]= max(abs(zall5_3Amm(i,:))); % find peak (positive or negative value) in that 1s time window
    zall3_5Apeak (i)= zall5_3Amm(i,X); % extract this peak value and put it in a new matrix
end
zall3_5Apeak=zall3_5Apeak';

%
% %TRIGGER 3  - animal 6
zall3_6 = zeros(size(dF_trials3_6));
for i = 1:size(dF_trials3_6,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_6(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_6(i,ind)); % baseline period stdev
    zall3_6(i,:)=(dF_trials3_6(i,:) - zb)/zsd; % Z score per bin
    zall3_6A (i,:)= zall3_6 (i,E:F); % extract 1s time window
    zall3_6Amm (i,:)= movmean(zall3_6A(i,:), 10); % transform data into moving average of 10ms
    [~,X]= max(abs(zall3_6Amm(i,:))); % find peak (positive or negative value) in that 1s time window
    zall3_6Apeak (i)= zall3_6Amm(i,X); % extract this peak value and put it in a new matrix
end
zall3_6Apeak=zall3_6Apeak';


%% %% %%  Quantify changes as min/max z score value for bout offset
% per 1s bin
C=1; % beginning of baseline  (1 means first data point of the 5s baseline)
D=bltimewin*fs; % end of baseline
% time unit is in ms. Multiply value by fs to obtain seconds.
% This example use a baseline of 5 seconds prior to trigger

BASELINE_Z = [C D]; 

E=bltimewin*fs; 
F=(bltimewin+1)*fs;

%TRIGGER 4  - animal 1
zall4_1 = zeros(size(dF_trials4_1));
for i = 1:size(dF_trials4_1,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials4_1(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials4_1(i,ind)); % baseline period stdev
    zall4_1(i,:)=(dF_trials4_1(i,:) - zb)/zsd; % Z score per bin
    zall4_1A (i,:)= zall4_1 (i,E:F); % extract 1s time window
    zall4_1Amm (i,:)= movmean(zall4_1A(i,:), 10); % transform data into moving average of 10ms
    [~,X]= max(abs(zall4_1Amm(i,:))); % find peak (positive or negative value) in that 1s time window
    zall4_1Apeak (i)= zall4_1Amm(i,X); % extract this peak value and put it in a new matrix
end
zall4_1Apeak=zall4_1Apeak';


%TRIGGER 4  - animal 2
zall4_2 = zeros(size(dF_trials4_2));
for i = 1:size(dF_trials4_2,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials4_2(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials4_2(i,ind)); % baseline period stdev
    zall4_2(i,:)=(dF_trials4_2(i,:) - zb)/zsd; % Z score per bin
    zall4_2A (i,:)= zall4_2 (i,E:F); % extract 1s time window
    zall4_2Amm (i,:)= movmean(zall4_2A(i,:), 10); % transform data into moving average of 10ms
    [~,X]= max(abs(zall4_2Amm(i,:))); % find peak (positive or negative value) in that 1s time window
    zall4_2Apeak (i)= zall4_2Amm(i,X); % extract this peak value and put it in a new matrix
end
zall4_2Apeak=zall4_2Apeak';


%TRIGGER 4  - animal 3
zall4_3 = zeros(size(dF_trials4_3));
for i = 1:size(dF_trials4_3,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials4_3(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials4_3(i,ind)); % baseline period stdev
    zall4_3(i,:)=(dF_trials4_3(i,:) - zb)/zsd; % Z score per bin
    zall4_3A (i,:)= zall4_3 (i,E:F); % extract 1s time window
    zall4_3Amm (i,:)= movmean(zall4_3A(i,:), 10); % transform data into moving average of 10ms
    [~,X]= max(abs(zall4_3Amm(i,:))); % find peak (positive or negative value) in that 1s time window
    zall4_3Apeak (i)= zall4_3Amm(i,X); % extract this peak value and put it in a new matrix
end
zall4_3Apeak=zall4_3Apeak';


%TRIGGER 4  - animal 4
zall4_4 = zeros(size(dF_trials4_4));
for i = 1:size(dF_trials4_4,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials4_4(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials4_4(i,ind)); % baseline period stdev
    zall4_4(i,:)=(dF_trials4_4(i,:) - zb)/zsd; % Z score per bin
    zall4_4A (i,:)= zall4_4 (i,E:F); % extract 1s time window
    zall4_4Amm (i,:)= movmean(zall4_4A(i,:), 10); % transform data into moving average of 10ms
    [~,X]= max(abs(zall4_4Amm(i,:))); % find peak (positive or negative value) in that 1s time window
    zall4_4Apeak (i)= zall4_4Amm(i,X); % extract this peak value and put it in a new matrix
end
zall4_4Apeak=zall4_4Apeak';


% %TRIGGER 4  - animal 5
zall4_5 = zeros(size(dF_trials4_5));
for i = 1:size(dF_trials4_5,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials4_5(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials4_5(i,ind)); % baseline period stdev
    zall4_5(i,:)=(dF_trials4_5(i,:) - zb)/zsd; % Z score per bin
    zall4_5A (i,:)= zall4_5(i,E:F); % extract 1s time window
    zall4_5Amm (i,:)= movmean(zall4_5A(i,:), 10); % transform data into moving average of 10ms
    [~,X]= max(abs(zall4_5Amm(i,:))); % find peak (positive or negative value) in that 1s time window
    zall4_5Apeak (i)= zall4_5Amm(i,X); % extract this peak value and put it in a new matrix
end
zall4_5Apeak=zall4_5Apeak';

%
% %TRIGGER 4  - animal 6
zall4_6 = zeros(size(dF_trials4_6));
for i = 1:size(dF_trials4_6,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials4_6(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials4_6(i,ind)); % baseline period stdev
    zall4_6(i,:)=(dF_trials4_6(i,:) - zb)/zsd; % Z score per bin
    zall4_6A (i,:)= zall4_6 (i,E:F); % extract 1s time window
    zall4_6Amm (i,:)= movmean(zall4_6A(i,:), 10); % transform data into moving average of 10ms
    [~,X]= max(abs(zall4_6Amm(i,:))); % find peak (positive or negative value) in that 1s time window
    zall4_6Apeak (i)= zall4_6Amm(i,X); % extract this peak value and put it in a new matrix
end
zall4_6Apeak=zall4_6Apeak';

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
