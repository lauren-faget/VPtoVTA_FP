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
file = 'VP GABA - VTA Glut - grp8';
experiment = 'Passive stim';
folder = 'D:\Hnasko lab_Lauren Faget_Data\Hnasko lab_Fiber Photometry\RESUME_FP_VP ChR2-VTA GCAMP\z-score Fig 300dpi';

%input notes for saved data
note_condition1 = '5Hz - 5s';
note_condition2 = '10Hz - 5s';
note_condition3 = '20Hz - 5s';
note_condition4 = '40Hz - 5s';

%% MODIFY trigger names here
%trigger names can be found in data/ epocs/ 
%format is '3 letters - underscore' or 'x - 3 numbers/numbers-letters combo - underscore'

EPOC_1 = 'x05H_'; % enter name of trigger 1
EPOC_2 = 'x10H_'; % enter name of trigger 2
EPOC_3 = 'x20H_'; % enter name of trigger 3
EPOC_4 = 'x40H_'; % enter name of trigger 4

%% Import data - Enter tank file paths here

fp1='tankfilepath1';
data1 = TDTbin2mat(fp1);

fp2='tankfilepath2';
data2 = TDTbin2mat(fp2);

fp3='tankfilepath3';
data3 = TDTbin2mat(fp3);

fp4='tankfilepath4';
data4 = TDTbin2mat(fp4);

fp5='tankfilepath5';
data5 = TDTbin2mat(fp5);

fp6='tankfilepath6';
data6 = TDTbin2mat(fp6);

fp7='tankfilepath7';
data7 = TDTbin2mat(fp7);

fp8='tankfilepath8';
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

%%  this will plot the non-normalized Gcamp/AF signals together
% figure;
% hold on;
% plot(Gc1,'k');
% plot(af1,'b');
% plot(cf1,'r');

%% B. Trigger timepoint extraction 
% plot the entire dF and the triggers you made

tr1_1=data1.epocs.(EPOC_1).onset;
tr2_1=data1.epocs.(EPOC_2).onset;
tr3_1=data1.epocs.(EPOC_3).onset;
tr4_1=data1.epocs.(EPOC_4).onset;

tr1_2=data2.epocs.(EPOC_1).onset;
tr2_2=data2.epocs.(EPOC_2).onset;
tr3_2=data2.epocs.(EPOC_3).onset;
tr4_2=data2.epocs.(EPOC_4).onset;

tr1_3=data3.epocs.(EPOC_1).onset;
tr2_3=data3.epocs.(EPOC_2).onset;
tr3_3=data3.epocs.(EPOC_3).onset;
tr4_3=data3.epocs.(EPOC_4).onset;

tr1_4=data4.epocs.(EPOC_1).onset;
tr2_4=data4.epocs.(EPOC_2).onset;
tr3_4=data4.epocs.(EPOC_3).onset;
tr4_4=data4.epocs.(EPOC_4).onset;

tr1_5=data5.epocs.(EPOC_1).onset;
tr2_5=data5.epocs.(EPOC_2).onset;
tr3_5=data5.epocs.(EPOC_3).onset;
tr4_5=data5.epocs.(EPOC_4).onset;

tr1_6=data6.epocs.(EPOC_1).onset;
tr2_6=data6.epocs.(EPOC_2).onset;
tr3_6=data6.epocs.(EPOC_3).onset;
tr4_6=data6.epocs.(EPOC_4).onset;

tr1_7=data7.epocs.(EPOC_1).onset;
tr2_7=data7.epocs.(EPOC_2).onset;
tr3_7=data7.epocs.(EPOC_3).onset;
tr4_7=data7.epocs.(EPOC_4).onset;

tr1_8=data8.epocs.(EPOC_1).onset;
tr2_8=data8.epocs.(EPOC_2).onset;
tr3_8=data8.epocs.(EPOC_3).onset;
tr4_8=data8.epocs.(EPOC_4).onset;

% 
%% B.this plots the entire dF and the triggers you made.

trt1_1=(tr1_1*fs);
trt2_1=(tr2_1*fs);
trt3_1=(tr3_1*fs);
trt4_1=(tr4_1*fs);

trt1_2=(tr1_2*fs);
trt2_2=(tr2_2*fs);
trt3_2=(tr3_2*fs);
trt4_2=(tr4_2*fs);

trt1_3=(tr1_3*fs);
trt2_3=(tr2_3*fs);
trt3_3=(tr3_3*fs);
trt4_3=(tr4_3*fs);

trt1_4=(tr1_4*fs);
trt2_4=(tr2_4*fs);
trt3_4=(tr3_4*fs);
trt4_4=(tr4_4*fs);

trt1_5=(tr1_5*fs);
trt2_5=(tr2_5*fs);
trt3_5=(tr3_5*fs);
trt4_5=(tr4_5*fs);

trt1_6=(tr1_6*fs);
trt2_6=(tr2_6*fs);
trt3_6=(tr3_6*fs);
trt4_6=(tr4_6*fs);

trt1_7=(tr1_7*fs);
trt2_7=(tr2_7*fs);
trt3_7=(tr3_7*fs);
trt4_7=(tr4_7*fs);

trt1_8=(tr1_8*fs);
trt2_8=(tr2_8*fs);
trt3_8=(tr3_8*fs);
trt4_8=(tr4_8*fs);


figure;
title('Entire Experiment - animal 1');hold on;
plot((1/fs:1/fs:length(dF1)/fs),dF1,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr1_1)
     plot([tr1_1(i1);tr1_1(i1)],[-15;15],'c');
 end
 for i2 = 1:length(tr2_1)
     plot([tr2_1(i2);tr2_1(i2)],[-15;15],'g');
 end
  for i3 = 1:length(tr3_1)
     plot([tr3_1(i3);tr3_1(i3)],[-15;15],'b');
 end
 for i4 = 1:length(tr4_1)
     plot([tr4_1(i4);tr4_1(i4)],[-15;15],'r');
 end

 figure;
title('Entire Experiment - animal 2');hold on;
plot((1/fs:1/fs:length(dF2)/fs),dF2,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

 for i1 = 1:length(tr1_2)
     plot([tr1_2(i1);tr1_2(i1)],[-15;15],'c');
 end
 for i2 = 1:length(tr2_2)
     plot([tr2_2(i2);tr2_2(i2)],[-15;15],'g');
 end
 for i3 = 1:length(tr3_2)
     plot([tr3_2(i3);tr3_2(i3)],[-15;15],'b');
 end
 for i4 = 1:length(tr4_2)
     plot([tr4_2(i4);tr4_2(i4)],[-15;15],'r');
 end

figure;
title('Entire Experiment - animal 3');hold on;
plot((1/fs:1/fs:length(dF3)/fs),dF3,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr1_3)
     plot([tr1_3(i1);tr1_3(i1)],[-15;15],'c');
 end
 for i2 = 1:length(tr2_3)
     plot([tr2_3(i2);tr2_3(i2)],[-15;15],'g');
 end
  for i3 = 1:length(tr3_3)
     plot([tr3_3(i3);tr3_3(i3)],[-15;15],'b');
 end
 for i4 = 1:length(tr4_3)
     plot([tr4_3(i4);tr4_3(i4)],[-15;15],'r');
 end
 
figure;
title('Entire Experiment - animal 4');hold on;
plot((1/fs:1/fs:length(dF4)/fs),dF4,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr1_4)
     plot([tr1_4(i1);tr1_4(i1)],[-15;15],'c');
 end
 for i2 = 1:length(tr2_4)
     plot([tr2_4(i2);tr2_4(i2)],[-15;15],'g');
 end
  for i3 = 1:length(tr3_4)
     plot([tr3_4(i3);tr3_4(i3)],[-15;15],'b');
 end
 for i4 = 1:length(tr4_4)
     plot([tr4_4(i4);tr4_4(i4)],[-15;15],'r');
 end

figure;
title('Entire Experiment - animal 5');hold on;
plot((1/fs:1/fs:length(dF5)/fs),dF5,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr1_5)
     plot([tr1_5(i1);tr1_5(i1)],[-15;15],'c');
 end
 for i2 = 1:length(tr2_5)
     plot([tr2_5(i2);tr2_5(i2)],[-15;15],'g');
 end
  for i3 = 1:length(tr3_5)
     plot([tr3_5(i3);tr3_5(i3)],[-15;15],'b');
 end
 for i4 = 1:length(tr4_5)
     plot([tr4_5(i4);tr4_5(i4)],[-15;15],'r');
 end

 figure;
title('Entire Experiment - animal 6');hold on;
plot((1/fs:1/fs:length(dF6)/fs),dF6,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr1_6)
     plot([tr1_6(i1);tr1_6(i1)],[-15;15],'c');
 end
 for i2 = 1:length(tr2_6)
     plot([tr2_6(i2);tr2_6(i2)],[-15;15],'g');
 end
  for i3 = 1:length(tr3_6)
     plot([tr3_6(i3);tr3_6(i3)],[-15;15],'b');
 end
 for i4 = 1:length(tr4_6)
     plot([tr4_6(i4);tr4_6(i4)],[-15;15],'r');
 end

  figure;
title('Entire Experiment - animal 7');hold on;
plot((1/fs:1/fs:length(dF7)/fs),dF7,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr1_7)
     plot([tr1_7(i1);tr1_7(i1)],[-15;15],'c');
 end
 for i2 = 1:length(tr2_7)
     plot([tr2_7(i2);tr2_7(i2)],[-15;15],'g');
 end
  for i3 = 1:length(tr3_7)
     plot([tr3_7(i3);tr3_7(i3)],[-15;15],'b');
 end
 for i4 = 1:length(tr4_7)
     plot([tr4_7(i4);tr4_7(i4)],[-15;15],'r');
 end

    figure;
title('Entire Experiment - animal 8');hold on;
plot((1/fs:1/fs:length(dF8)/fs),dF8,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr1_8)
     plot([tr1_8(i1);tr1_8(i1)],[-15;15],'c');
 end
 for i2 = 1:length(tr2_8)
     plot([tr2_8(i2);tr2_8(i2)],[-15;15],'g');
 end
  for i3 = 1:length(tr3_8)
     plot([tr3_8(i3);tr3_8(i3)],[-15;15],'b');
 end
 for i4 = 1:length(tr4_8)
     plot([tr4_8(i4);tr4_8(i4)],[-15;15],'r');
 end


% close all;

%% C.this plots +/- time windows around the triggers made during experiment.

%DEFINE your time window here;
bltimewin=5;  %% baseline time window in seconds
timewin=10;   %% time window after trigger in seconds
timescale = bltimewin + timewin;

dF_trials1_1=[]; 
dF_trials2_1=[]; 
dF_trials3_1=[]; 
dF_trials4_1=[]; 

dF_trials1_2=[]; 
dF_trials2_2=[]; 
dF_trials3_2=[]; 
dF_trials4_2=[]; 

dF_trials1_3=[]; 
dF_trials2_3=[]; 
dF_trials3_3=[]; 
dF_trials4_3=[]; 

dF_trials1_4=[]; 
dF_trials2_4=[]; 
dF_trials3_4=[]; 
dF_trials4_4=[];

dF_trials1_5=[]; 
dF_trials2_5=[]; 
dF_trials3_5=[]; 
dF_trials4_5=[]; 

dF_trials1_6=[]; 
dF_trials2_6=[]; 
dF_trials3_6=[]; 
dF_trials4_6=[]; 

dF_trials1_7=[]; 
dF_trials2_7=[]; 
dF_trials3_7=[]; 
dF_trials4_7=[]; 

dF_trials1_8=[]; 
dF_trials2_8=[]; 
dF_trials3_8=[]; 
dF_trials4_8=[]; 


%animal 1
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt1_1)
     if (trt1_1(i1)+timewin*fs<length(dF1) && trt1_1(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials1_1(count,:)=dF1(trt1_1(i1)-(bltimewin*fs):trt1_1(i1)+(timewin*fs));  
    end
end
% extract trials for the second trigger.
count=0;
for i1 = 1:length(trt2_1)
     if (trt2_1(i1)+timewin*fs<length(dF1) && trt2_1(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials2_1(count,:)=dF1(trt2_1(i1)-(bltimewin*fs):trt2_1(i1)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i3 = 1:length(trt3_1)
    if (trt3_1(i3)+timewin*fs<length(dF1) && trt3_1(i3)-timewin*fs>0)     
        count=count+1;
        dF_trials3_1(count,:)=dF1(trt3_1(i3)-(bltimewin*fs):trt3_1(i3)+(timewin*fs));       
    end
end

% extract trials for the fourth trigger.
count=0;
for i4 = 1:length(trt4_1)
    if (trt4_1(i4)+timewin*fs<length(dF1) && trt4_1(i4)-timewin*fs>0)     
        count=count+1;
        dF_trials4_1(count,:)=dF1(trt4_1(i4)-(bltimewin*fs):trt4_1(i4)+(timewin*fs));       
    end
end

%animal 2
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt1_2)
     if (trt1_2(i1)+timewin*fs<length(dF2) && trt1_2(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials1_2(count,:)=dF2(trt1_2(i1)-(bltimewin*fs):trt1_2(i1)+(timewin*fs));  
    end
end
% extract trials for the second trigger.
count=0;
for i1 = 1:length(trt2_2)
     if (trt2_2(i1)+timewin*fs<length(dF2) && trt2_2(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials2_2(count,:)=dF2(trt2_2(i1)-(bltimewin*fs):trt2_2(i1)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i3 = 1:length(trt3_2)
    if (trt3_2(i3)+timewin*fs<length(dF2) && trt3_2(i3)-timewin*fs>0)     
        count=count+1;
        dF_trials3_2(count,:)=dF2(trt3_2(i3)-(bltimewin*fs):trt3_2(i3)+(timewin*fs));       
    end
end
% extract trials for the fourth trigger.
count=0;
for i4 = 1:length(trt4_2)
    if (trt4_2(i4)+timewin*fs<length(dF2) && trt4_2(i4)-timewin*fs>0)     
        count=count+1;
        dF_trials4_2(count,:)=dF2(trt4_2(i4)-(bltimewin*fs):trt4_2(i4)+(timewin*fs));       
    end
end

%animal 3
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt1_3)
     if (trt1_3(i1)+timewin*fs<length(dF3) && trt1_3(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials1_3(count,:)=dF3(trt1_3(i1)-(bltimewin*fs):trt1_3(i1)+(timewin*fs));  
    end
end
% extract trials for the second trigger.
count=0;
for i1 = 1:length(trt2_3)
     if (trt2_3(i1)+timewin*fs<length(dF3) && trt2_3(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials2_3(count,:)=dF3(trt2_3(i1)-(bltimewin*fs):trt2_3(i1)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i3 = 1:length(trt3_3)
    if (trt3_3(i3)+timewin*fs<length(dF3) && trt3_3(i3)-timewin*fs>0)     
        count=count+1;
        dF_trials3_3(count,:)=dF3(trt3_3(i3)-(bltimewin*fs):trt3_3(i3)+(timewin*fs));       
    end
end
% extract trials for the fourth trigger.
count=0;
for i4 = 1:length(trt4_3)
    if (trt4_3(i4)+timewin*fs<length(dF3) && trt4_3(i4)-timewin*fs>0)     
        count=count+1;
        dF_trials4_3(count,:)=dF3(trt4_3(i4)-(bltimewin*fs):trt4_3(i4)+(timewin*fs));       
    end
end

%animal 4
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt1_4)
     if (trt1_4(i1)+timewin*fs<length(dF4) && trt1_4(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials1_4(count,:)=dF4(trt1_4(i1)-(bltimewin*fs):trt1_4(i1)+(timewin*fs));  
    end
end
% extract trials for the second trigger.
count=0;
for i1 = 1:length(trt2_4)
     if (trt2_4(i1)+timewin*fs<length(dF4) && trt2_4(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials2_4(count,:)=dF4(trt2_4(i1)-(bltimewin*fs):trt2_4(i1)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i3 = 1:length(trt3_4)
    if (trt3_4(i3)+timewin*fs<length(dF4) && trt3_4(i3)-timewin*fs>0)     
        count=count+1;
        dF_trials3_4(count,:)=dF4(trt3_4(i3)-(bltimewin*fs):trt3_4(i3)+(timewin*fs));       
    end
end
% extract trials for the fourth trigger.
count=0;
for i4 = 1:length(trt4_4)
    if (trt4_4(i4)+timewin*fs<length(dF4) && trt4_4(i4)-timewin*fs>0)     
        count=count+1;
        dF_trials4_4(count,:)=dF4(trt4_4(i4)-(bltimewin*fs):trt4_4(i4)+(timewin*fs));       
    end
end

%animal 5
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt1_5)
     if (trt1_5(i1)+timewin*fs<length(dF5) && trt1_5(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials1_5(count,:)=dF5(trt1_5(i1)-(bltimewin*fs):trt1_5(i1)+(timewin*fs));  
    end
end
% extract trials for the second trigger.
count=0;
for i1 = 1:length(trt2_5)
     if (trt2_5(i1)+timewin*fs<length(dF5) && trt2_5(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials2_5(count,:)=dF5(trt2_5(i1)-(bltimewin*fs):trt2_5(i1)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i3 = 1:length(trt3_5)
    if (trt3_5(i3)+timewin*fs<length(dF5) && trt3_5(i3)-timewin*fs>0)     
        count=count+1;
        dF_trials3_5(count,:)=dF5(trt3_5(i3)-(bltimewin*fs):trt3_5(i3)+(timewin*fs));       
    end
end
% extract trials for the fourth trigger.
count=0;
for i4 = 1:length(trt4_5)
    if (trt4_5(i4)+timewin*fs<length(dF5) && trt4_5(i4)-timewin*fs>0)     
        count=count+1;
        dF_trials4_5(count,:)=dF5(trt4_5(i4)-(bltimewin*fs):trt4_5(i4)+(timewin*fs));       
    end
end

%animal 6
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt1_6)
     if (trt1_6(i1)+timewin*fs<length(dF6) && trt1_6(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials1_6(count,:)=dF6(trt1_6(i1)-(bltimewin*fs):trt1_6(i1)+(timewin*fs));  
    end
end
% extract trials for the second trigger.
count=0;
for i1 = 1:length(trt2_6)
     if (trt2_6(i1)+timewin*fs<length(dF6) && trt2_6(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials2_6(count,:)=dF6(trt2_6(i1)-(bltimewin*fs):trt2_6(i1)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i3 = 1:length(trt3_6)
    if (trt3_6(i3)+timewin*fs<length(dF6) && trt3_6(i3)-timewin*fs>0)     
        count=count+1;
        dF_trials3_6(count,:)=dF6(trt3_6(i3)-(bltimewin*fs):trt3_6(i3)+(timewin*fs));       
    end
end
% extract trials for the fourth trigger.
count=0;
for i4 = 1:length(trt4_6)
    if (trt4_6(i4)+timewin*fs<length(dF6) && trt4_6(i4)-timewin*fs>0)     
        count=count+1;
        dF_trials4_6(count,:)=dF6(trt4_6(i4)-(bltimewin*fs):trt4_6(i4)+(timewin*fs));       
    end
end

%animal 7
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt1_7)
     if (trt1_7(i1)+timewin*fs<length(dF7) && trt1_7(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials1_7(count,:)=dF7(trt1_7(i1)-(bltimewin*fs):trt1_7(i1)+(timewin*fs));  
    end
end
% extract trials for the second trigger.
count=0;
for i1 = 1:length(trt2_7)
     if (trt2_7(i1)+timewin*fs<length(dF7) && trt2_7(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials2_7(count,:)=dF7(trt2_7(i1)-(bltimewin*fs):trt2_7(i1)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i3 = 1:length(trt3_7)
    if (trt3_7(i3)+timewin*fs<length(dF7) && trt3_7(i3)-timewin*fs>0)     
        count=count+1;
        dF_trials3_7(count,:)=dF7(trt3_7(i3)-(bltimewin*fs):trt3_7(i3)+(timewin*fs));       
    end
end
% extract trials for the fourth trigger.
count=0;
for i4 = 1:length(trt4_7)
    if (trt4_7(i4)+timewin*fs<length(dF7) && trt4_7(i4)-timewin*fs>0)     
        count=count+1;
        dF_trials4_7(count,:)=dF7(trt4_7(i4)-(bltimewin*fs):trt4_7(i4)+(timewin*fs));       
    end
end

%animal 8
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt1_8)
     if (trt1_8(i1)+timewin*fs<length(dF8) && trt1_8(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials1_8(count,:)=dF8(trt1_8(i1)-(bltimewin*fs):trt1_8(i1)+(timewin*fs));  
    end
end
% extract trials for the second trigger.
count=0;
for i1 = 1:length(trt2_8)
     if (trt2_8(i1)+timewin*fs<length(dF8) && trt2_8(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials2_8(count,:)=dF8(trt2_8(i1)-(bltimewin*fs):trt2_8(i1)+(timewin*fs));  
    end
end
% extract trials for the third trigger.
count=0;
for i3 = 1:length(trt3_8)
    if (trt3_8(i3)+timewin*fs<length(dF8) && trt3_8(i3)-timewin*fs>0)     
        count=count+1;
        dF_trials3_8(count,:)=dF8(trt3_8(i3)-(bltimewin*fs):trt3_8(i3)+(timewin*fs));       
    end
end

% extract trials for the fourth trigger.
count=0;
for i4 = 1:length(trt4_8)
    if (trt4_8(i4)+timewin*fs<length(dF8) && trt4_8(i4)-timewin*fs>0)     
        count=count+1;
        dF_trials4_8(count,:)=dF8(trt4_8(i4)-(bltimewin*fs):trt4_8(i4)+(timewin*fs));       
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

bl2_1=dF_trials2_1(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm2_1=mean(bl2_1,2);
dFn2_1=(dF_trials2_1-repmat(blm2_1,1,size(dF_trials2_1,2)));

bl3_1=dF_trials3_1(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_1=mean(bl3_1,2);
dFn3_1=(dF_trials3_1-repmat(blm3_1,1,size(dF_trials3_1,2)));

bl4_1=dF_trials4_1(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm4_1=mean(bl4_1,2);
dFn4_1=(dF_trials4_1-repmat(blm4_1,1,size(dF_trials4_1,2)));


bl1_2=dF_trials1_2(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_2=mean(bl1_2,2);
dFn1_2=(dF_trials1_2-repmat(blm1_2,1,size(dF_trials1_2,2)));

bl2_2=dF_trials2_2(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm2_2=mean(bl2_2,2);
dFn2_2=(dF_trials2_2-repmat(blm2_2,1,size(dF_trials2_2,2)));

bl3_2=dF_trials3_2(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_2=mean(bl3_2,2);
dFn3_2=(dF_trials3_2-repmat(blm3_2,1,size(dF_trials3_2,2)));

bl4_2=dF_trials4_2(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm4_2=mean(bl4_2,2);
dFn4_2=(dF_trials4_2-repmat(blm4_2,1,size(dF_trials4_2,2)));



bl1_3=dF_trials1_3(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_3=mean(bl1_3,2);
dFn1_3=(dF_trials1_3-repmat(blm1_3,1,size(dF_trials1_3,2)));

bl2_3=dF_trials2_3(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm2_3=mean(bl2_3,2);
dFn2_3=(dF_trials2_3-repmat(blm2_3,1,size(dF_trials2_3,2)));

bl3_3=dF_trials3_3(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_3=mean(bl3_3,2);
dFn3_3=(dF_trials3_3-repmat(blm3_3,1,size(dF_trials3_3,2)));

bl4_3=dF_trials4_3(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm4_3=mean(bl4_3,2);
dFn4_3=(dF_trials4_3-repmat(blm4_3,1,size(dF_trials4_3,2)));



bl1_4=dF_trials1_4(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_4=mean(bl1_4,2);
dFn1_4=(dF_trials1_4-repmat(blm1_4,1,size(dF_trials1_4,2)));

bl2_4=dF_trials2_4(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm2_4=mean(bl2_4,2);
dFn2_4=(dF_trials2_4-repmat(blm2_4,1,size(dF_trials2_4,2)));

bl3_4=dF_trials3_4(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_4=mean(bl3_4,2);
dFn3_4=(dF_trials3_4-repmat(blm3_4,1,size(dF_trials3_4,2)));

bl4_4=dF_trials4_4(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm4_4=mean(bl4_4,2);
dFn4_4=(dF_trials4_4-repmat(blm4_4,1,size(dF_trials4_4,2)));


bl1_5=dF_trials1_5(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_5=mean(bl1_5,2);
dFn1_5=(dF_trials1_5-repmat(blm1_5,1,size(dF_trials1_5,2)));

bl2_5=dF_trials2_5(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm2_5=mean(bl2_5,2);
dFn2_5=(dF_trials2_5-repmat(blm2_5,1,size(dF_trials2_5,2)));

bl3_5=dF_trials3_5(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_5=mean(bl3_5,2);
dFn3_5=(dF_trials3_5-repmat(blm3_5,1,size(dF_trials3_5,2)));

bl4_5=dF_trials4_5(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm4_5=mean(bl4_5,2);
dFn4_5=(dF_trials4_5-repmat(blm4_5,1,size(dF_trials4_5,2)));



bl1_6=dF_trials1_6(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_6=mean(bl1_6,2);
dFn1_6=(dF_trials1_6-repmat(blm1_6,1,size(dF_trials1_6,2)));

bl2_6=dF_trials2_6(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm2_6=mean(bl2_6,2);
dFn2_6=(dF_trials2_6-repmat(blm2_6,1,size(dF_trials2_6,2)));

bl3_6=dF_trials3_6(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_6=mean(bl3_6,2);
dFn3_6=(dF_trials3_6-repmat(blm3_6,1,size(dF_trials3_6,2)));

bl4_6=dF_trials4_6(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm4_6=mean(bl4_6,2);
dFn4_6=(dF_trials4_6-repmat(blm4_6,1,size(dF_trials4_6,2)));


bl1_7=dF_trials1_7(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_7=mean(bl1_7,2);
dFn1_7=(dF_trials1_7-repmat(blm1_7,1,size(dF_trials1_7,2)));

bl2_7=dF_trials2_7(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm2_7=mean(bl2_7,2);
dFn2_7=(dF_trials2_7-repmat(blm2_7,1,size(dF_trials2_7,2)));

bl3_7=dF_trials3_7(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_7=mean(bl3_7,2);
dFn3_7=(dF_trials3_7-repmat(blm3_7,1,size(dF_trials3_7,2)));

bl4_7=dF_trials4_7(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm4_7=mean(bl4_7,2);
dFn4_7=(dF_trials4_7-repmat(blm4_7,1,size(dF_trials4_7,2)));


bl1_8=dF_trials1_8(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm1_8=mean(bl1_8,2);
dFn1_8=(dF_trials1_8-repmat(blm1_8,1,size(dF_trials1_8,2)));

bl2_8=dF_trials2_8(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm2_8=mean(bl2_8,2);
dFn2_8=(dF_trials2_8-repmat(blm2_8,1,size(dF_trials2_8,2)));

bl3_8=dF_trials3_8(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_8=mean(bl3_8,2);
dFn3_8=(dF_trials3_8-repmat(blm3_8,1,size(dF_trials3_8,2)));

bl4_8=dF_trials4_8(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm4_8=mean(bl4_8,2);
dFn4_8=(dF_trials4_8-repmat(blm4_8,1,size(dF_trials4_8,2)));


%Non-normalized dF/F for all trials and all animals
dF_trials1= [dF_trials1_1;dF_trials1_2;dF_trials1_3;dF_trials1_4;dF_trials1_5;dF_trials1_6;dF_trials1_7;dF_trials1_8];
dF_trials2= [dF_trials2_1;dF_trials2_2;dF_trials2_3;dF_trials2_4;dF_trials2_5;dF_trials2_6;dF_trials2_7;dF_trials2_8];
dF_trials3= [dF_trials3_1;dF_trials3_2;dF_trials3_3;dF_trials3_4;dF_trials3_5;dF_trials3_6;dF_trials3_7;dF_trials3_8];
dF_trials4= [dF_trials4_1;dF_trials4_2;dF_trials4_3;dF_trials4_4;dF_trials4_5;dF_trials4_6;dF_trials4_7;dF_trials4_8];

%Normalized dF/F for all trials and all animals
dFn1= [dFn1_1;dFn1_2;dFn1_3;dFn1_4;dFn1_5;dFn1_6;dFn1_7;dFn1_8];
dFn2= [dFn2_1;dFn2_2;dFn2_3;dFn2_4;dFn2_5;dFn2_6;dFn2_7;dFn2_8];
dFn3= [dFn3_1;dFn3_2;dFn3_3;dFn3_4;dFn3_5;dFn3_6;dFn3_7;dFn3_8];
dFn4= [dFn4_1;dFn4_2;dFn4_3;dFn4_4;dFn4_5;dFn4_6;dFn4_7;dFn4_8];


figure;
imagesc(dFn1); colorbar; hold on;
caxis([-5 20]);  % Y axis scale
xline(5*fs,'-k','linewidth',2); % create a line at trigger onset
xticks((0:fs:timescale*fs));  % X axis scale
xticklabels({-5 : 10});  % X axis scale labeling
title('normalized dF/F trigger 1');
ylabel('Trials', 'FontSize', 12);
xlabel('Time (s)','FontSize',12)
set (gca,'TickDir','out');

figure;
imagesc(dFn2); colorbar; hold on;
caxis([-5 20]);
xline(5*fs,'-k','linewidth',2);
xticks((0:fs:timescale*fs));
xticklabels({-5 : 10}) ;
title('normalized dF/F trigger 2');
ylabel('Trials', 'FontSize', 12);
xlabel('Time (s)','FontSize',12)
set (gca,'TickDir','out');

figure;
imagesc(dFn3); colorbar; hold on;
caxis([-5 20]);
xline(5*fs,'-k','linewidth',2);
xticks((0:fs:timescale*fs));
xticklabels({-5 : 10});
title('normalized dF/F trigger 3');
ylabel('Trials', 'FontSize', 12);
xlabel('Time (s)','FontSize',12)
set (gca,'TickDir','out');

figure;
imagesc(dFn4); colorbar; hold on;
caxis([-5 20]);
xline(5*fs,'-k','linewidth',2);
xticks((0:fs:timescale*fs));
xticklabels({-5 : 10});
title('normalized dF/F trigger 4');
ylabel('Trials', 'FontSize', 12);
xlabel('Time (s)','FontSize',12)
set (gca,'TickDir','out');

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
% subplot(223)
% imagesc(dFn2); colorbar; title('Normalized data trigger 2');
% caxis([-5 20]);
% xticks((0:fs:timescale*fs));
% xticklabels({-5 : 10}); hold on;
% subplot(224);
% imagesc(dF_trials2);colorbar;title('Non-Normalized data trigger 2');
% caxis([-5 20]);
% xticks((0:fs:timescale*fs));
% xticklabels({-5 : 10});


%% Normalized Mean +/- SEM with shadded error boars

dFmnnm1=mean(dFn1(1:end,1:end));
dFsemnm1=std(dFn1(1:end,1:end))/sqrt(size(dFn1,1));
dFmnnm2=mean(dFn2(1:end,1:end));
dFsemnm2=std(dFn2(1:end,1:end))/sqrt(size(dFn2,1));
dFmnnm3=mean(dFn3(1:end,1:end));
dFsemnm3=std(dFn3(1:end,1:end))/sqrt(size(dFn3,1));
dFmnnm4=mean(dFn4(1:end,1:end));
dFsemnm4=std(dFn4(1:end,1:end))/sqrt(size(dFn4,1));

figure;
shadedErrorBar((-bltimewin:1/fs:timewin),dFmnnm1,dFsemnm1); 
title('normalized dF/F trigger 1');
axis([-bltimewin inf -5 20]);
hold on;
plot([0;0],[-5;20],'b','Linewidth',1);
plot([1;1],[-5;20],':b','Linewidth',1); % plot a second line at the end of the stim
ylabel('dF/F', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');

figure;
shadedErrorBar((-bltimewin:1/fs:timewin),dFmnnm2,dFsemnm2); 
title('normalized dF/F trigger 2');
axis([-bltimewin inf -5 20]);
hold on;
plot([0;0],[-5;20],'b','Linewidth',1);
plot([1;1],[-5;20],':b','Linewidth',1);
ylabel('dF/F', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');

figure;
shadedErrorBar((-bltimewin:1/fs:timewin),dFmnnm3,dFsemnm3); 
title('normalized dF/F trigger 3');
axis([-bltimewin inf -5 20]);
hold on;
plot([0;0],[-5;20],'b','Linewidth',1);
plot([1;1],[-5;20],':b','Linewidth',1);
ylabel('dF/F', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');

figure;
shadedErrorBar((-bltimewin:1/fs:timewin),dFmnnm4,dFsemnm4); 
title('normalized dF/F trigger 4');
axis([-bltimewin inf -5 20]);
hold on;
plot([0;0],[-5;20],'b','Linewidth',1);
plot([1;1],[-5;20],':b','Linewidth',1);
ylabel('dF/F', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');

%% plotting triggers against each other normalized
figure;
hold on;
plot((-bltimewin:1/fs:timewin),dFmnnm1,'c');
plot((-bltimewin:1/fs:timewin),dFmnnm2,'g');
plot((-bltimewin:1/fs:timewin),dFmnnm3,'r');
plot((-bltimewin:1/fs:timewin),dFmnnm4,'b');
axis([-bltimewin inf -5 20]);
hold on;
plot([0;0],[-5;20],'k','Linewidth',1);
plot([1;1],[-5;20],':k','Linewidth',1);
legend('trigger 1','trigger 2','trigger 3','trigger 4');
title('Normalized DF/F laser');
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

%TRIGGER 2
zall2 = zeros(size(dF_trials2));
for i = 1:size(dF_trials2,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials2(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials2(i,ind)); % baseline period stdev
    zall2(i,:)=(dF_trials2(i,:) - zb)/zsd; % Z score per bin
end
zerror2 = std(zall2)/sqrt(size(zall2,1));
zmean2 =  mean (zall2);

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

%TRIGGER 4
zall4 = zeros(size(dF_trials4));
for i = 1:size(dF_trials4,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials4(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials4(i,ind)); % baseline period stdev
    zall4(i,:)=(dF_trials4(i,:) - zb)/zsd; % Z score per bin
end
zerror4 = std(zall4)/sqrt(size(zall4,1));
zmean4 =  mean (zall4);

%FIGURES
% heat map - colorbar graph - trigger 1
figure; 
imagesc(zall1); colorbar; hold on;
caxis([-2 5]);
xticks((0:fs:timescale*fs));
xticklabels({-5 : 10});
ylabel('Trials', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition1,'-','zscore colorbar'); 
title(note_title);

% save figure as matlab fig and pdf
cd (folder);
saveas(gcf,note_title, 'fig');
saveas(gcf,note_title, 'pdf');

% mean +/- sem shadded error bar graph - trigger 1
figure; 
shadedErrorBar([-bltimewin:1/fs:timewin],zmean1,zerror1, 'lineProps','b'); hold on;
axis([-5 inf -2 5]);
plot([0;0],[-2;5],'k','Linewidth',1);
ylabel('z score', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition1,'-','zscore'); 
title(note_title);

% save figure as matlab fig and pdf
cd (folder);
saveas(gcf,note_title, 'fig');
saveas(gcf,note_title, 'pdf');

% heat map - colorbar graph - trigger 2
figure; 
imagesc(zall2); colorbar; hold on;
caxis([-2 5]);
xticks((0:fs:timescale*fs));
xticklabels({-5 : 10});
title('Z-Score Trigger 2 Heat Map')
ylabel('Trials', 'FontSize', 12);
xlabel('Time (s)','FontSize',12); 
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition2,'-','zscore colorbar'); 
title(note_title);

% save figure as matlab fig and pdf
cd (folder);
saveas(gcf,note_title, 'fig');
saveas(gcf,note_title, 'pdf');

% mean +/- sem shadded error bar graph - trigger 2
figure; 
shadedErrorBar([-bltimewin:1/fs:timewin],zmean2,zerror2, 'lineProps','b'); hold on;
axis([-5 inf -2 5]);
plot([0;0],[-2;5],'k','Linewidth',1);
ylabel('z score', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition2,'-','zscore'); 
title(note_title);

% save figure as matlab fig and pdf
cd (folder);
saveas(gcf,note_title, 'fig');
saveas(gcf,note_title, 'pdf');



% heat map - colorbar graph - trigger 3
figure;
imagesc(zall3); colorbar; hold on;
caxis([-2 5]);
xticks((0:fs:timescale*fs));
xticklabels({-5 : 10});
ylabel('Trials', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition3,'-','zscore colorbar'); 
title(note_title);

% save figure as matlab fig and pdf
cd (folder);
saveas(gcf,note_title, 'fig');
saveas(gcf,note_title, 'pdf');

% mean +/- sem shadded error bar graph - trigger 3
figure; 
shadedErrorBar([-bltimewin:1/fs:timewin],zmean3,zerror3, 'lineProps','b'); hold on;
axis([-5 inf -2 5]);
plot([0;0],[-2;5],'k','Linewidth',1);
ylabel('z score', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition3,'-','zscore'); 
title(note_title);

% save figure as matlab fig and pdf
cd (folder);
saveas(gcf,note_title, 'fig');
saveas(gcf,note_title, 'pdf');



% heat map - colorbar graph - trigger 4
figure;
imagesc(zall4); colorbar; hold on;
caxis([-2 5]);
xticks((0:fs:timescale*fs));
xticklabels({-5 : 10});
ylabel('Trials', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition4,'-','zscore colorbar'); 
title(note_title);

% save figure as matlab fig and pdf
cd (folder);
saveas(gcf,note_title, 'fig');
saveas(gcf,note_title, 'pdf');

% mean +/- sem shadded error bar graph - trigger 4
figure; 
shadedErrorBar([-bltimewin:1/fs:timewin],zmean4,zerror4, 'lineProps','b'); hold on;
axis([-5 inf -2 5]);
plot([0;0],[-2;5],'k','Linewidth',1);
ylabel('z score', 'FontSize', 12);
xlabel('Time (s)','FontSize',12);
set (gca,'TickDir','out');
note_title =  strcat(file,'-',experiment,'-',note_condition4,'-','zscore'); 
title(note_title);

% save figure as matlab fig and pdf
cd (folder);
saveas(gcf,note_title, 'fig');
saveas(gcf,note_title, 'pdf');


%%  Quantify changes as area under the curve - AUC

%Modify time window as needed - example with 1s prior stim and 1s during

% TRIGGER 1
AUC1=[]; % 
AUC1(1,1)=trapz(zmean1(:,(bltimewin-1)*fs:(bltimewin)*fs));
AUC1(1,2)=trapz(zmean1(:,(bltimewin)*fs:(bltimewin+1)*fs));

figure;
hBar = bar(AUC1, 'FaceColor', [.8 .8 .8]);
ylim([-500 2500]); hold on;
% centers = get(hBar, 'XData');
% plot(centers(1:2), [1 1]*AUC(1,2)*1.1, '-k', 'LineWidth', 2)
% p1 = plot(mean(centers(1:2)), AUC(1,2)*1.2, '*k');
set(gca,'xticklabel',{'baseline','stim'});
set (gca,'TickDir','out');
title({ 'AUC trigger 1'});


% TRIGGER 2
AUC2=[]; % 
AUC2(1,1)=trapz(zmean2(:,(bltimewin-1)*fs:(bltimewin)*fs));
AUC2(1,2)=trapz(zmean2(:,(bltimewin)*fs:(bltimewin+1)*fs));

figure;
hBar = bar(AUC2, 'FaceColor', [.8 .8 .8]);
ylim([-500 2500]); hold on;
% centers = get(hBar, 'XData');
% plot(centers(1:2), [1 1]*AUC(1,2)*1.1, '-k', 'LineWidth', 2)
% p1 = plot(mean(centers(1:2)), AUC(1,2)*1.2, '*k');
set(gca,'xticklabel',{'baseline','stim'});
set (gca,'TickDir','out');
title({ 'AUC trigger 2'})

% TRIGGER 3
AUC3=[]; % 
AUC3(1,1)=trapz(zmean3(:,(bltimewin-1)*fs:(bltimewin)*fs));
AUC3(1,2)=trapz(zmean3(:,(bltimewin)*fs:(bltimewin+1)*fs));

figure;
hBar = bar(AUC3, 'FaceColor', [.8 .8 .8]);
ylim([-500 2500]); hold on;
% centers = get(hBar, 'XData');
% plot(centers(1:2), [1 1]*AUC(1,2)*1.1, '-k', 'LineWidth', 2)
% p1 = plot(mean(centers(1:2)), AUC(1,2)*1.2, '*k');
set(gca,'xticklabel',{'baseline','stim'});
set (gca,'TickDir','out');
title({ 'AUC trigger 3'})

% TRIGGER 4
AUC4=[]; % 
AUC4(1,1)=trapz(zmean4(:,(bltimewin-1)*fs:(bltimewin)*fs));
AUC4(1,2)=trapz(zmean4(:,(bltimewin)*fs:(bltimewin+1)*fs));

figure;
hBar = bar(AUC4, 'FaceColor', [.8 .8 .8]);
ylim([-500 2500]); hold on;
% centers = get(hBar, 'XData');
% plot(centers(1:2), [1 1]*AUC(1,2)*1.1, '-k', 'LineWidth', 2)
% p1 = plot(mean(centers(1:2)), AUC(1,2)*1.2, '*k');
set(gca,'xticklabel',{'baseline','stim'});
set (gca,'TickDir','out');
title({ 'AUC trigger 4'})


%% Save matlab data file
% You can save space by saving only necessary data - deleting unwanted data in the
% workspace - before saving .mat file

% cd 'filepath';
% save filename.mat
