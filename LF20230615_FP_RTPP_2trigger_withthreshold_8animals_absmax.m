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
experiment = 'RTPP - day 1';
folder = 'filepath';

%input notes for saved data
note_condition = 'paired onset 2s off - 5s on';

%% MODIFY trigger names here
%trigger names can be found in data/ epocs/ 
%format is '3 letters - underscore' or 'x - 3 numbers/numbers-letters combo - underscore'

EPOC_1 = 'pai_'; % enter name of trigger 1


%% DEFINE time threshold in paired and unpaired side around trigger here

% minimum time after trigger
TIME_THRESHOLD = 5; %  in seconds 

% minimum time before trigger - BASELINE
TIME_THRESHOLD2 = 2; %  in seconds 


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


%% B. Trigger timepoint extraction - onset & offset
% plot the entire dF and the triggers you made

tr1_1=data1.epocs.(EPOC_1).onset;
tr2_1=data1.epocs.(EPOC_1).offset;

tr1_2=data2.epocs.(EPOC_1).onset;
tr2_2=data2.epocs.(EPOC_1).offset;

tr1_3=data3.epocs.(EPOC_1).onset;
tr2_3=data3.epocs.(EPOC_1).offset;

tr1_4=data4.epocs.(EPOC_1).onset;
tr2_4=data4.epocs.(EPOC_1).offset;

tr1_5=data5.epocs.(EPOC_1).onset;
tr2_5=data5.epocs.(EPOC_1).offset;

tr1_6=data6.epocs.(EPOC_1).onset;
tr2_6=data6.epocs.(EPOC_1).offset;

tr1_7=data7.epocs.(EPOC_1).onset;
tr2_7=data7.epocs.(EPOC_1).offset;

tr1_8=data8.epocs.(EPOC_1).onset;
tr2_8=data8.epocs.(EPOC_1).offset;

% 
%% %% %% Define new array of long periods in paired side
% laser on - data 1

laser_on1 = (tr1_1);
laser_off1 = [0 ; tr2_1];
laser_on_duration1 = laser_off1 (2:end) - laser_on1 ;
laser_off_duration1 = laser_on1 (1:end) - laser_off1 (1:end-1) ;

average_laser_on_duration1 = mean(laser_on_duration1);
average_laser_off_duration1 = mean(laser_off_duration1);

laser_on_duration_ind1 = find(laser_on_duration1 >= TIME_THRESHOLD);
laser_off_duration_ind1 = find(laser_off_duration1 >= TIME_THRESHOLD2);
laser_common_ind1 = intersect (laser_on_duration_ind1,laser_off_duration_ind1);
% 

diff_ind_i1 = 1;
for i1 = 1:length(laser_common_ind1)
    data1.epocs.las_.onset(i1) = laser_on1(diff_ind_i1);
    data1.epocs.las_.data(i1) = 1; % set the data value, arbitrary 1
    diff_ind_i1 = laser_common_ind1(i1); % increment the index
end

% Transpose the arrays to make them column vectors like other epocs
data1.epocs.las_.onset = laser_on1(laser_common_ind1)';
tr3_1=data1.epocs.las_.onset';

%% laser off - data 1

laser_on1 = (tr1_1);
laser_off1 = (tr2_1);
laser_on_duration1 = laser_off1 - laser_on1 ;
laser_off_duration1 = laser_on1 (2:end) - laser_off1 (1:end-1) ;
    
laser_on_duration_ind1 = find(laser_on_duration1 >= TIME_THRESHOLD2);
laser_off_duration_ind1 = find(laser_off_duration1 >= TIME_THRESHOLD);
laser_common_ind1 = intersect (laser_on_duration_ind1,laser_off_duration_ind1);

% 
diff_off_ind_i2 = 1;
for i2 = 1:length(laser_common_ind1)
    data1.epocs.lasoff_.onset(i2) = laser_off1(diff_off_ind_i2);
    data1.epocs.lasoff_.data(i2) = 1; % set the data value, arbitrary 1
    diff_off_ind_i2 = laser_common_ind1(i2); % increment the index
end


% Transpose the arrays to make them column vectors like other epocs
data1.epocs.lasoff_.onset = laser_off1(laser_common_ind1)';
tr4_1=data1.epocs.lasoff_.onset';

%% %% Define new array of long periods in paired side - data 2
% laser on - data 2

laser_on2 = (tr1_2);
laser_off2 = [0 ; tr2_2];
laser_on_duration2 = laser_off2 (2:end) - laser_on2 ;
laser_off_duration2 = laser_on2 (1:end) - laser_off2 (1:end-1) ;

average_laser_on_duration2 = mean(laser_on_duration2);
average_laser_off_duration2 = mean(laser_off_duration2);

laser_on_duration_ind2 = find(laser_on_duration2 >= TIME_THRESHOLD);
laser_off_duration_ind2 = find(laser_off_duration2 >= TIME_THRESHOLD2);
laser_common_ind2 = intersect (laser_on_duration_ind2,laser_off_duration_ind2);

% 
diff_ind_i1 = 1;
for i1 = 1:length(laser_common_ind2)
    data2.epocs.las_.onset(i1) = laser_on2(diff_ind_i1);
    data2.epocs.las_.data(i1) = 1; % set the data value, arbitrary 1
    diff_ind_i1 = laser_common_ind2(i1); % increment the index
end

% Transpose the arrays to make them column vectors like other epocs
data2.epocs.las_.onset = laser_on2(laser_common_ind2)';
tr3_2=data2.epocs.las_.onset';

%% laser off - data 2
 
laser_on2 = (tr1_2);
laser_off2 = (tr2_2);
laser_on_duration2 = laser_off2 - laser_on2 ;
laser_off_duration2 = laser_on2 (2:end) - laser_off2 (1:end-1) ;
    
laser_on_duration_ind2 = find(laser_on_duration2 >= TIME_THRESHOLD2);
laser_off_duration_ind2 = find(laser_off_duration2 >= TIME_THRESHOLD);
laser_common_ind2 = intersect (laser_on_duration_ind2,laser_off_duration_ind2);
% 
diff_off_ind_i2 = 1;
for i2 = 1:length(laser_common_ind2)
    data2.epocs.lasoff_.onset(i2) = laser_off2(diff_off_ind_i2);
    data2.epocs.lasoff_.data(i2) = 1; % set the data value, arbitrary 1
    diff_off_ind_i2 = laser_common_ind2(i2); % increment the index
end


% Transpose the arrays to make them column vectors like other epocs
data2.epocs.lasoff_.onset = laser_off2(laser_common_ind2)';
tr4_2=data2.epocs.lasoff_.onset';

%% %% Define new array of long periods in paired side - data 3
% laser on - data 3

laser_on3 = (tr1_3);
laser_off3 = [0 ; tr2_3];
laser_on_duration3 = laser_off3 (2:end) - laser_on3 ;
laser_off_duration3 = laser_on3 (1:end) - laser_off3 (1:end-1) ;

average_laser_on_duration3 = mean(laser_on_duration3);
average_laser_off_duration3 = mean(laser_off_duration3);

laser_on_duration_ind3 = find(laser_on_duration3 >= TIME_THRESHOLD);
laser_off_duration_ind3 = find(laser_off_duration3 >= TIME_THRESHOLD2);
laser_common_ind3 = intersect (laser_on_duration_ind3,laser_off_duration_ind3);
% 
diff_ind_i1 = 1;
for i1 = 1:length(laser_common_ind3)
    data3.epocs.las_.onset(i1) = laser_on3(diff_ind_i1);
    data3.epocs.las_.data(i1) = 1; % set the data value, arbitrary 1
    diff_ind_i1 = laser_common_ind3(i1); % increment the index
end

% Transpose the arrays to make them column vectors like other epocs
data3.epocs.las_.onset = laser_on3(laser_common_ind3)';
tr3_3=data3.epocs.las_.onset';

%% laser off - data 3
    
laser_on3 = (tr1_3);
laser_off3 = (tr2_3);
laser_on_duration3 = laser_off3 - laser_on3 ;
laser_off_duration3 = laser_on3 (2:end) - laser_off3 (1:end-1) ;
    
laser_on_duration_ind3 = find(laser_on_duration3 >= TIME_THRESHOLD2);
laser_off_duration_ind3 = find(laser_off_duration3 >= TIME_THRESHOLD);
laser_common_ind3 = intersect (laser_on_duration_ind3,laser_off_duration_ind3);
% 
diff_off_ind_i2 = 1;
for i2 = 1:length(laser_common_ind3)
    data3.epocs.lasoff_.onset(i2) = laser_off3(diff_off_ind_i2);
    data3.epocs.lasoff_.data(i2) = 1; % set the data value, arbitrary 1
    diff_off_ind_i2 = laser_common_ind3(i2); % increment the index
end


% Transpose the arrays to make them column vectors like other epocs
data3.epocs.lasoff_.onset = laser_off3(laser_common_ind3)';
tr4_3=data3.epocs.lasoff_.onset';

%% %% Define new array of long periods in paired side - data 4
% laser on - data 4

laser_on4 = (tr1_4);
laser_off4 = [0 ; tr2_4];
laser_on_duration4 = laser_off4 (2:end) - laser_on4 ;
laser_off_duration4 = laser_on4 (1:end) - laser_off4 (1:end-1) ;

average_laser_on_duration4 = mean(laser_on_duration4);
average_laser_off_duration4 = mean(laser_off_duration4);

laser_on_duration_ind4 = find(laser_on_duration4 >= TIME_THRESHOLD);
laser_off_duration_ind4 = find(laser_off_duration4 >= TIME_THRESHOLD2);
laser_common_ind4 = intersect (laser_on_duration_ind4,laser_off_duration_ind4);
% 
diff_ind_i1 = 1;
for i1 = 1:length(laser_common_ind4)
    data4.epocs.las_.onset(i1) = laser_on4(diff_ind_i1);
    data4.epocs.las_.data(i1) = 1; % set the data value, arbitrary 1
    diff_ind_i1 = laser_common_ind4(i1); % increment the index
end

% Transpose the arrays to make them column vectors like other epocs
data4.epocs.las_.onset = laser_on4(laser_common_ind4)';
tr3_4=data4.epocs.las_.onset';

%% laser off - data 4
    
laser_on4 = (tr1_4);
laser_off4 = (tr2_4);
laser_on_duration4 = laser_off4 - laser_on4 ;
laser_off_duration4 = laser_on4 (2:end) - laser_off4 (1:end-1) ;
    
laser_on_duration_ind4 = find(laser_on_duration4 >= TIME_THRESHOLD2);
laser_off_duration_ind4 = find(laser_off_duration4 >= TIME_THRESHOLD);
laser_common_ind4 = intersect (laser_on_duration_ind4,laser_off_duration_ind4);

% 
diff_off_ind_i2 = 1;
for i2 = 1:length(laser_common_ind4)
    data4.epocs.lasoff_.onset(i2) = laser_off4(diff_off_ind_i2);
    data4.epocs.lasoff_.data(i2) = 1; % set the data value, arbitrary 1
    diff_off_ind_i2 = laser_common_ind4(i2); % increment the index
end


% Transpose the arrays to make them column vectors like other epocs
data4.epocs.lasoff_.onset = laser_off4(laser_common_ind4)';
tr4_4=data4.epocs.lasoff_.onset';

%% %% Define new array of long periods in paired side - data 5
% laser on - data 5

laser_on5 = (tr1_5);
laser_off5 = [0 ; tr2_5];
laser_on_duration5 = laser_off5 (2:end) - laser_on5 ;
laser_off_duration5 = laser_on5 (1:end) - laser_off5 (1:end-1) ;

average_laser_on_duration5 = mean(laser_on_duration5);
average_laser_off_duration5 = mean(laser_off_duration5);

laser_on_duration_ind5 = find(laser_on_duration5 >= TIME_THRESHOLD);
laser_off_duration_ind5 = find(laser_off_duration5 >= TIME_THRESHOLD2);
laser_common_ind5 = intersect (laser_on_duration_ind5,laser_off_duration_ind5);
% 
diff_ind_i1 = 1;
for i1 = 1:length(laser_common_ind5)
    data5.epocs.las_.onset(i1) = laser_on5(diff_ind_i1);
    data5.epocs.las_.data(i1) = 1; % set the data value, arbitrary 1
    diff_ind_i1 = laser_common_ind5(i1); % increment the index
end

% Transpose the arrays to make them column vectors like other epocs
data5.epocs.las_.onset = laser_on5(laser_common_ind5)';
tr3_5=data5.epocs.las_.onset';

%% laser off - data 5
    
laser_on5 = (tr1_5);
laser_off5 = (tr2_5);
laser_on_duration5 = laser_off5 - laser_on5 ;
laser_off_duration5 = laser_on5 (2:end) - laser_off5 (1:end-1) ;
    
laser_on_duration_ind5 = find(laser_on_duration5 >= TIME_THRESHOLD2);
laser_off_duration_ind5 = find(laser_off_duration5 >= TIME_THRESHOLD);
laser_common_ind5 = intersect (laser_on_duration_ind5,laser_off_duration_ind5);

% 
diff_off_ind_i2 = 1;
for i2 = 1:length(laser_common_ind5)
    data5.epocs.lasoff_.onset(i2) = laser_off5(diff_off_ind_i2);
    data5.epocs.lasoff_.data(i2) = 1; % set the data value, arbitrary 1
    diff_off_ind_i2 = laser_common_ind5(i2); % increment the index
end


% Transpose the arrays to make them column vectors like other epocs
data5.epocs.lasoff_.onset = laser_off5(laser_common_ind5)';
tr4_5=data5.epocs.lasoff_.onset';

%% %% Define new array of long periods in paired side - data 6
% laser on - data 6

laser_on6 = (tr1_6);
laser_off6 = [0 ; tr2_6];
laser_on_duration6 = laser_off6 (2:end) - laser_on6 ;
laser_off_duration6 = laser_on6 (1:end) - laser_off6 (1:end-1) ;

average_laser_on_duration6 = mean(laser_on_duration6);
average_laser_off_duration6 = mean(laser_off_duration6);

laser_on_duration_ind6 = find(laser_on_duration6 >= TIME_THRESHOLD);
laser_off_duration_ind6 = find(laser_off_duration6 >= TIME_THRESHOLD2);
laser_common_ind6 = intersect (laser_on_duration_ind6,laser_off_duration_ind6);
% 
diff_ind_i1 = 1;
for i1 = 1:length(laser_common_ind6)
    data6.epocs.las_.onset(i1) = laser_on6(diff_ind_i1);
    data6.epocs.las_.data(i1) = 1; % set the data value, arbitrary 1
    diff_ind_i1 = laser_common_ind6(i1); % increment the index
end

% Transpose the arrays to make them column vectors like other epocs
data6.epocs.las_.onset = laser_on6(laser_common_ind6)';
tr3_6=data6.epocs.las_.onset';

%% laser off - data 6
    
laser_on6 = (tr1_6);
laser_off6 = (tr2_6);
laser_on_duration6 = laser_off6 - laser_on6 ;
laser_off_duration6 = laser_on6 (2:end) - laser_off6 (1:end-1) ;
    
laser_on_duration_ind6 = find(laser_on_duration6 >= TIME_THRESHOLD2);
laser_off_duration_ind6 = find(laser_off_duration6 >= TIME_THRESHOLD);
laser_common_ind6 = intersect (laser_on_duration_ind6,laser_off_duration_ind6);
% 
diff_off_ind_i2 = 1;
for i2 = 1:length(laser_common_ind6)
    data6.epocs.lasoff_.onset(i2) = laser_off6(diff_off_ind_i2);
    data6.epocs.lasoff_.data(i2) = 1; % set the data value, arbitrary 1
    diff_off_ind_i2 = laser_common_ind6(i2); % increment the index
end


% Transpose the arrays to make them column vectors like other epocs
data6.epocs.lasoff_.onset = laser_off6(laser_common_ind6)';
tr4_6=data6.epocs.lasoff_.onset';

%% %% Define new array of long periods in paired side - data 7
% laser on - data 7

laser_on7 = (tr1_7);
laser_off7 = [0 ; tr2_7];
laser_on_duration7 = laser_off7 (2:end) - laser_on7 ;
laser_off_duration7 = laser_on7 (1:end) - laser_off7 (1:end-1) ;

average_laser_on_duration7 = mean(laser_on_duration7);
average_laser_off_duration7 = mean(laser_off_duration7);

laser_on_duration_ind7 = find(laser_on_duration7 >= TIME_THRESHOLD);
laser_off_duration_ind7 = find(laser_off_duration7 >= TIME_THRESHOLD2);
laser_common_ind7 = intersect (laser_on_duration_ind7,laser_off_duration_ind7);
% 
diff_ind_i1 = 1;
for i1 = 1:length(laser_common_ind7)
    data7.epocs.las_.onset(i1) = laser_on7(diff_ind_i1);
    data7.epocs.las_.data(i1) = 1; % set the data value, arbitrary 1
    diff_ind_i1 = laser_common_ind7(i1); % increment the index
end

% Transpose the arrays to make them column vectors like other epocs
data7.epocs.las_.onset = laser_on7(laser_common_ind7)';
tr3_7=data7.epocs.las_.onset';

%% laser off - data 7
    
laser_on7 = (tr1_7);
laser_off7 = (tr2_7);
laser_on_duration7 = laser_off7 - laser_on7 ;
laser_off_duration7 = laser_on7 (2:end) - laser_off7 (1:end-1) ;
    
laser_on_duration_ind7 = find(laser_on_duration7 >= TIME_THRESHOLD2);
laser_off_duration_ind7 = find(laser_off_duration7 >= TIME_THRESHOLD);
laser_common_ind7 = intersect (laser_on_duration_ind7,laser_off_duration_ind7);
% 
diff_off_ind_i2 = 1;
for i2 = 1:length(laser_common_ind7)
    data7.epocs.lasoff_.onset(i2) = laser_off7(diff_off_ind_i2);
    data7.epocs.lasoff_.data(i2) = 1; % set the data value, arbitrary 1
    diff_off_ind_i2 = laser_common_ind7(i2); % increment the index
end


% Transpose the arrays to make them column vectors like other epocs
data7.epocs.lasoff_.onset = laser_off7(laser_common_ind7)';
tr4_7=data7.epocs.lasoff_.onset';

%% %% Define new array of long periods in paired side - data 8
% laser on - data 8

laser_on8 = (tr1_8);
laser_off8 = [0 ; tr2_8];
laser_on_duration8 = laser_off8 (2:end) - laser_on8 ;
laser_off_duration8 = laser_on8 (1:end) - laser_off8 (1:end-1) ;

average_laser_on_duration8 = mean(laser_on_duration8);
average_laser_off_duration8 = mean(laser_off_duration8);

laser_on_duration_ind8 = find(laser_on_duration8 >= TIME_THRESHOLD);
laser_off_duration_ind8 = find(laser_off_duration8 >= TIME_THRESHOLD2);
laser_common_ind8 = intersect (laser_on_duration_ind8,laser_off_duration_ind8);
% 
diff_ind_i1 = 1;
for i1 = 1:length(laser_common_ind8)
    data8.epocs.las_.onset(i1) = laser_on8(diff_ind_i1);
    data8.epocs.las_.data(i1) = 1; % set the data value, arbitrary 1
    diff_ind_i1 = laser_common_ind8(i1); % increment the index
end

% Transpose the arrays to make them column vectors like other epocs
data8.epocs.las_.onset = laser_on8(laser_common_ind8)';
tr3_8=data8.epocs.las_.onset';

%% laser off - data 8
    
laser_on8 = (tr1_8);
laser_off8 = (tr2_8);
laser_on_duration8 = laser_off8 - laser_on8 ;
laser_off_duration8 = laser_on8 (2:end) - laser_off8 (1:end-1) ;
    
laser_on_duration_ind8 = find(laser_on_duration8 >= TIME_THRESHOLD2);
laser_off_duration_ind8 = find(laser_off_duration8 >= TIME_THRESHOLD);
laser_common_ind8 = intersect (laser_on_duration_ind8,laser_off_duration_ind8);
% 
diff_off_ind_i2 = 1;
for i2 = 1:length(laser_common_ind8)
    data8.epocs.lasoff_.onset(i2) = laser_off8(diff_off_ind_i2);
    data8.epocs.lasoff_.data(i2) = 1; % set the data value, arbitrary 1
    diff_off_ind_i2 = laser_common_ind8(i2); % increment the index
end


% Transpose the arrays to make them column vectors like other epocs
data8.epocs.lasoff_.onset = laser_off8(laser_common_ind8)';
tr4_8=data8.epocs.lasoff_.onset';

%% B.this plots the entire dF and the triggers you made.

trt3_1=(tr3_1*fs);

trt3_2=(tr3_2*fs);

trt3_3=(tr3_3*fs);

trt3_4=(tr3_4*fs);

trt3_5=(tr3_5*fs);

trt3_6=(tr3_6*fs);

trt3_7=(tr3_7*fs);

trt3_8=(tr3_8*fs);


figure;
title('Entire Experiment - animal 1');hold on;
plot((1/fs:1/fs:length(dF1)/fs),dF1,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr3_1)
     plot([tr3_1(i1);tr3_1(i1)],[-15;15],'c');
 end

 figure;
title('Entire Experiment - animal 2');hold on;
plot((1/fs:1/fs:length(dF2)/fs),dF2,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

 for i1 = 1:length(tr3_2)
     plot([tr3_2(i1);tr3_2(i1)],[-15;15],'c');
 end

figure;
title('Entire Experiment - animal 3');hold on;
plot((1/fs:1/fs:length(dF3)/fs),dF3,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr3_3)
     plot([tr3_3(i1);tr3_3(i1)],[-15;15],'c');
 end
 
figure;
title('Entire Experiment - animal 4');hold on;
plot((1/fs:1/fs:length(dF4)/fs),dF4,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr3_4)
     plot([tr3_4(i1);tr3_4(i1)],[-15;15],'c');
 end

figure;
title('Entire Experiment - animal 5');hold on;
plot((1/fs:1/fs:length(dF5)/fs),dF5,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr3_5)
     plot([tr3_5(i1);tr3_5(i1)],[-15;15],'c');
 end

 figure;
title('Entire Experiment - animal 6');hold on;
plot((1/fs:1/fs:length(dF6)/fs),dF6,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr3_6)
     plot([tr3_6(i1);tr3_6(i1)],[-15;15],'c');
 end

  figure;
title('Entire Experiment - animal 7');hold on;
plot((1/fs:1/fs:length(dF7)/fs),dF7,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr3_7)
     plot([tr3_7(i1);tr3_7(i1)],[-15;15],'c');
 end

    figure;
title('Entire Experiment - animal 8');hold on;
plot((1/fs:1/fs:length(dF8)/fs),dF8,'k'); hold on;
axis ([-inf inf -15 15]); 
hold on;

  for i1 = 1:length(tr3_8)
     plot([tr3_8(i1);tr3_8(i1)],[-15;15],'c');
 end

close all;

%% C.this plots +/- time windows around the triggers made during experiment.

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

dF_trials3_7=[]; 

dF_trials3_8=[]; 


%animal 1
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt3_1)
     if (trt3_1(i1)+timewin*fs<length(dF1) && trt3_1(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials3_1(count,:)=dF1(trt3_1(i1)-(bltimewin*fs):trt3_1(i1)+(timewin*fs));  
    end
end

%animal 2
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt3_2)
     if (trt3_2(i1)+timewin*fs<length(dF2) && trt3_2(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials3_2(count,:)=dF2(trt3_2(i1)-(bltimewin*fs):trt3_2(i1)+(timewin*fs));  
    end
end

%animal 3
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt3_3)
     if (trt3_3(i1)+timewin*fs<length(dF3) && trt3_3(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials3_3(count,:)=dF3(trt3_3(i1)-(bltimewin*fs):trt3_3(i1)+(timewin*fs));  
    end
end

%animal 4
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt3_4)
     if (trt3_4(i1)+timewin*fs<length(dF4) && trt3_4(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials3_4(count,:)=dF4(trt3_4(i1)-(bltimewin*fs):trt3_4(i1)+(timewin*fs));  
    end
end

%animal 5
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt3_5)
     if (trt3_5(i1)+timewin*fs<length(dF5) && trt3_5(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials3_5(count,:)=dF5(trt3_5(i1)-(bltimewin*fs):trt3_5(i1)+(timewin*fs));  
    end
end

%animal 6
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt3_6)
     if (trt3_6(i1)+timewin*fs<length(dF6) && trt3_6(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials3_6(count,:)=dF6(trt3_6(i1)-(bltimewin*fs):trt3_6(i1)+(timewin*fs));  
    end
end

%animal 7
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt3_7)
     if (trt3_7(i1)+timewin*fs<length(dF7) && trt3_7(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials3_7(count,:)=dF7(trt3_7(i1)-(bltimewin*fs):trt3_7(i1)+(timewin*fs));  
    end
end

%animal 8
% extract trials for the first trigger.
count=0;
for i1 = 1:length(trt3_8)
     if (trt3_8(i1)+timewin*fs<length(dF8) && trt3_8(i1)-timewin*fs>0)     
        count=count+1;
        dF_trials3_8(count,:)=dF8(trt3_8(i1)-(bltimewin*fs):trt3_8(i1)+(timewin*fs));  
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

bl3_1=dF_trials3_1(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_1=mean(bl3_1,2);
dFn3_1=(dF_trials3_1-repmat(blm3_1,1,size(dF_trials3_1,2)));

bl3_2=dF_trials3_2(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_2=mean(bl3_2,2);
dFn3_2=(dF_trials3_2-repmat(blm3_2,1,size(dF_trials3_2,2)));

bl3_3=dF_trials3_3(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_3=mean(bl3_3,2);
dFn3_3=(dF_trials3_3-repmat(blm3_3,1,size(dF_trials3_3,2)));

bl3_4=dF_trials3_4(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_4=mean(bl3_4,2);
dFn3_4=(dF_trials3_4-repmat(blm3_4,1,size(dF_trials3_4,2)));

bl3_5=dF_trials3_5(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_5=mean(bl3_5,2);
dFn3_5=(dF_trials3_5-repmat(blm3_5,1,size(dF_trials3_5,2)));

bl3_6=dF_trials3_6(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_6=mean(bl3_6,2);
dFn3_6=(dF_trials3_6-repmat(blm3_6,1,size(dF_trials3_6,2)));

bl3_7=dF_trials3_7(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_7=mean(bl3_7,2);
dFn3_7=(dF_trials3_7-repmat(blm3_7,1,size(dF_trials3_7,2)));

bl3_8=dF_trials3_8(:,BASELINE_NORM (1):BASELINE_NORM (2)); 
blm3_8=mean(bl3_8,2);
dFn3_8=(dF_trials3_8-repmat(blm3_8,1,size(dF_trials3_8,2)));



%Non-normalized dF/F for all trials and all animals
dF_trials3= [dF_trials3_1;dF_trials3_2;dF_trials3_3;dF_trials3_4;dF_trials3_5;dF_trials3_6;dF_trials3_7;dF_trials3_8];

%Normalized dF/F for all trials and all animals
dFn3= [dFn3_1;dFn3_2;dFn3_3;dFn3_4;dFn3_5;dFn3_6;dFn3_7;dFn3_8];


figure;
imagesc(dFn3); colorbar; hold on;
caxis([-12 20]);  % Y axis scale
xline(5*fs,'-k','linewidth',2); % create a line at trigger onset
xticks((0:fs:timescale*fs));  % X axis scale
xticklabels({-5 : 10});  % X axis scale labeling
title('normalized dF/F trigger 3');
ylabel('Trials', 'FontSize', 12);
xlabel('Time (s)','FontSize',12)
set (gca,'TickDir','out');


%% Normalized Mean +/- SEM with shadded error boars

dFmnnm3=mean(dFn3(1:end,1:end));
dFsemnm3=std(dFn3(1:end,1:end))/sqrt(size(dFn3,1));

figure;
shadedErrorBar((-bltimewin:1/fs:timewin),dFmnnm3,dFsemnm3); 
title('normalized dF/F trigger 3');
axis([-bltimewin inf -12 20]);
hold on;
plot([0;0],[-12;20],'b','Linewidth',1);
% plot([1;1],[-12;20],':b','Linewidth',1); % plot a second line at the end of the stim
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
% heat map - colorbar graph - trigger 1
figure; 
imagesc(zall3); colorbar; hold on;
caxis([-12 20]);
xticks((0:fs:timescale*fs));
xticklabels({-5 : 10});
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
shadedErrorBar([-bltimewin:1/fs:timewin],zmean3,zerror3, 'lineProps','b'); hold on;
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
%%  Quantify changes as min/max z score value
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

%moving average per 10 ms
zmean3_1 = movmean(zmean3_1a, 10);

zmean3_1A= zmean3_1 (:,E:F);
zmean3_1B= zmean3_1 (:,F:G);
zmean3_1C= zmean3_1 (:,H:I);
zmean3_1D= zmean3_1 (:,I:J);
zmean3_1E= zmean3_1 (:,G:K);
zmean3_1F= zmean3_1 (:,L:E);
zmean3_1G= zmean3_1 (:,J:M);
zmean3_1H= zmean3_1 (:,K:N);

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

%moving average per 10 ms
zmean3_2 = movmean(zmean3_2a, 10);

zmean3_2A= zmean3_2 (:,E:F);
zmean3_2B= zmean3_2 (:,F:G);
zmean3_2C= zmean3_2 (:,H:I);
zmean3_2D= zmean3_2 (:,I:J);
zmean3_2E= zmean3_2 (:,G:K);
zmean3_2F= zmean3_2 (:,L:E);
zmean3_2G= zmean3_2 (:,J:M);
zmean3_2H= zmean3_2 (:,K:N);

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

%moving average per 10 ms
zmean3_3 = movmean(zmean3_3a, 10);

zmean3_3A= zmean3_3 (:,E:F);
zmean3_3B= zmean3_3 (:,F:G);
zmean3_3C= zmean3_3 (:,H:I);
zmean3_3D= zmean3_3 (:,I:J);
zmean3_3E= zmean3_3 (:,G:K);
zmean3_3F= zmean3_3 (:,L:E);
zmean3_3G= zmean3_3 (:,J:M);
zmean3_3H= zmean3_3 (:,K:N);

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

%moving average per 10 ms
zmean3_4 = movmean(zmean3_4a, 10);

zmean3_4A= zmean3_4 (:,E:F);
zmean3_4B= zmean3_4 (:,F:G);
zmean3_4C= zmean3_4 (:,H:I);
zmean3_4D= zmean3_4 (:,I:J);
zmean3_4E= zmean3_4 (:,G:K);
zmean3_4F= zmean3_4 (:,L:E);
zmean3_4G= zmean3_4 (:,J:M);
zmean3_4H= zmean3_4 (:,K:N);

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
% 
%moving average per 10 ms
zmean3_5 = movmean(zmean3_5a, 10);

zmean3_5A= zmean3_5 (:,E:F);
zmean3_5B= zmean3_5 (:,F:G);
zmean3_5C= zmean3_5 (:,H:I);
zmean3_5D= zmean3_5 (:,I:J);
zmean3_5E= zmean3_5 (:,G:K);
zmean3_5F= zmean3_5 (:,L:E);
zmean3_5G= zmean3_5 (:,J:M);
zmean3_5H= zmean3_5 (:,K:N);

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
%
%moving average per 10 ms
zmean3_6 = movmean(zmean3_6a, 10);

zmean3_6A= zmean3_6 (:,E:F);
zmean3_6B= zmean3_6 (:,F:G);
zmean3_6C= zmean3_6 (:,H:I);
zmean3_6D= zmean3_6 (:,I:J);
zmean3_6E= zmean3_6 (:,G:K);
zmean3_6F= zmean3_6 (:,L:E);
zmean3_6G= zmean3_6 (:,J:M);
zmean3_6H= zmean3_6 (:,K:N);

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


% TRIGGER 3  - animal 7
zall3_7 = zeros(size(dF_trials3_7));
for i = 1:size(dF_trials3_7,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_7(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_7(i,ind)); % baseline period stdev
    zall3_7(i,:)=(dF_trials3_7(i,:) - zb)/zsd; % Z score per bin
end
zerror3_7 = std(zall3_7)/sqrt(size(zall3_7,1));
zmean3_7a =  mean (zall3_7);

%moving average per 10 ms
zmean3_7 = movmean(zmean3_7a, 10);

zmean3_7A= zmean3_7 (:,E:F);
zmean3_7B= zmean3_7 (:,F:G);
zmean3_7C= zmean3_7 (:,H:I);
zmean3_7D= zmean3_7 (:,I:J);
zmean3_7E= zmean3_7 (:,G:K);
zmean3_7F= zmean3_7 (:,L:E);
zmean3_7G= zmean3_7 (:,J:M);
zmean3_7H= zmean3_7 (:,K:N);

[~,X] = max(abs(zmean3_7A));
zmean3_7Apeak = zmean3_7A(X);
[~,X] = max(abs(zmean3_7B));
zmean3_7Bpeak = zmean3_7B(X);
[~,X] = max(abs(zmean3_7C));
zmean3_7Cpeak = zmean3_7C(X);
[~,X] = max(abs(zmean3_7D));
zmean3_7Dpeak = zmean3_7D(X);
[~,X] = max(abs(zmean3_7E));
zmean3_7Epeak = zmean3_7E(X);
[~,X] = max(abs(zmean3_7F));
zmean3_7Fpeak = zmean3_7F(X);
[~,X] = max(abs(zmean3_7G));
zmean3_7Gpeak = zmean3_7G(X);
[~,X] = max(abs(zmean3_7H));
zmean3_7Hpeak = zmean3_7H(X);


%TRIGGER 3  - animal 8
zall3_8 = zeros(size(dF_trials3_8));
for i = 1:size(dF_trials3_8,1)
    ind = (BASELINE_Z(1) :  BASELINE_Z(2));
    zb = mean(dF_trials3_8(i,ind)); % baseline period mean (-5sec to 0sec)
    zsd = std(dF_trials3_8(i,ind)); % baseline period stdev
    zall3_8(i,:)=(dF_trials3_8(i,:) - zb)/zsd; % Z score per bin
end
zerror3_8 = std(zall3_8)/sqrt(size(zall3_8,1));
zmean3_8a =  mean (zall3_8);

%moving average per 10 ms
zmean3_8 = movmean(zmean3_8a, 10);

zmean3_8A= zmean3_8 (:,E:F);
zmean3_8B= zmean3_8 (:,F:G);
zmean3_8C= zmean3_8 (:,H:I);
zmean3_8D= zmean3_8 (:,I:J);
zmean3_8E= zmean3_8 (:,G:K);
zmean3_8F= zmean3_8 (:,L:E);
zmean3_8G= zmean3_8 (:,J:M);
zmean3_8H= zmean3_8 (:,K:N);

[~,X] = max(abs(zmean3_8A));
zmean3_8Apeak = zmean3_8A(X);
[~,X] = max(abs(zmean3_8B));
zmean3_8Bpeak = zmean3_8B(X);
[~,X] = max(abs(zmean3_8C));
zmean3_8Cpeak = zmean3_8C(X);
[~,X] = max(abs(zmean3_8D));
zmean3_8Dpeak = zmean3_8D(X);
[~,X] = max(abs(zmean3_8E));
zmean3_8Epeak = zmean3_8E(X);
[~,X] = max(abs(zmean3_8F));
zmean3_8Fpeak = zmean3_8F(X);
[~,X] = max(abs(zmean3_8G));
zmean3_8Gpeak = zmean3_8G(X);
[~,X] = max(abs(zmean3_8H));
zmean3_8Hpeak = zmean3_8H(X);



zabsmax_PreStim21= [zmean3_1Fpeak;zmean3_2Fpeak;zmean3_3Fpeak;zmean3_4Fpeak;zmean3_5Fpeak;zmean3_6Fpeak;zmean3_7Fpeak;zmean3_8Fpeak];
zabsmax_Prestim10= [zmean3_1Apeak;zmean3_2Apeak;zmean3_3Apeak;zmean3_4Apeak;zmean3_5Apeak;zmean3_6Apeak;zmean3_7Apeak;zmean3_8Apeak];
zabsmax_PostStim01= [zmean3_1Bpeak;zmean3_2Bpeak;zmean3_3Bpeak;zmean3_4Bpeak;zmean3_5Bpeak;zmean3_6Bpeak;zmean3_7Bpeak;zmean3_8Bpeak];
zabsmax_PostStim12= [zmean3_1Epeak;zmean3_2Epeak;zmean3_3Epeak;zmean3_4Epeak;zmean3_5Epeak;zmean3_6Epeak;zmean3_7Epeak;zmean3_8Epeak];
zabsmax_PostStim23= [zmean3_1Hpeak;zmean3_2Hpeak;zmean3_3Hpeak;zmean3_4Hpeak;zmean3_5Hpeak;zmean3_6Hpeak;zmean3_7Hpeak;zmean3_8Hpeak];
zabsmax_PostStim45= [zmean3_1Cpeak;zmean3_2Cpeak;zmean3_3Cpeak;zmean3_4Cpeak;zmean3_5Cpeak;zmean3_6Cpeak;zmean3_7Cpeak;zmean3_8Cpeak];
zabsmax_PostStim56= [zmean3_1Dpeak;zmean3_2Dpeak;zmean3_3Dpeak;zmean3_4Dpeak;zmean3_5Dpeak;zmean3_6Dpeak;zmean3_7Dpeak;zmean3_8Dpeak];
zabsmax_PostStim67= [zmean3_1Gpeak;zmean3_2Gpeak;zmean3_3Gpeak;zmean3_4Gpeak;zmean3_5Gpeak;zmean3_6Gpeak;zmean3_7Gpeak;zmean3_8Gpeak];

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


%% Save matlab data file
% You can save space by saving only necessary data - deleting unwanted data in the
% workspace - before saving .mat file

% cd 'filepath';
% save filename.mat
