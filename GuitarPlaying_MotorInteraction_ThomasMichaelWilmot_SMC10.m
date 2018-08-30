%% Guitar Playing and Motor Interaction
% Thomas Michael Wilmot
% Student of Sound and Music Computing MSc at Aalborg Univeristy
% Last updated: 03/08/18
%
% This script is used to analyse electromyoneurograms recorded when a
% subject performed the simplified motion components that combine to create
% the strumming motion used in guitar playing.
%
% For more information contact: twilmo14@student.aau.dk

clc;
clear all;
close all;

%% Go to data folder

% Change directory
cd('C:\Users\twilm\OneDrive\Documents\GitHub\SMC10_ThomasWilmot\SMC10ThomasWilmot\Delsysdata210318\TomSimpleMotionDelsys_210318');

%% Import Data
% if file below is not shown then:
%cd('C:\Users\Thomas\Documents\GitHub\SMC10_ThomasWilmot\SMC10ThomasWilmot');
% manually open script called 'Reformat_Delsys_data_Thomas_Wilmot' > select
% all > press F9. Then re-run the code in the current script
load('DelsysData_Thomas_Wilmot_SMC10');

%% Normalise motion data (max and min values are given by Delsys)

for idx_var=1:length(zworksp_idx);
    name_var=(zworksp_idx(idx_var).name);
    curr_var=eval(name_var);
    % Acc(+/-16)
    curr_var(idx_Acc,:)=curr_var(idx_Acc,:)/16;
    % Gyro(+/-2000)
    curr_var(idx_Gyro,:)=curr_var(idx_Gyro,:)/2000;
    % Mag(+/-1000)
    curr_var(idx_Mag,:)=curr_var(idx_Mag,:)/1000;
    assignin('base',(zworksp_idx(idx_var).name),curr_var);
end
clear('idx_var','curr_var','name_var');

%% Find indices for hand position data

% Indices for dorsum gyroscope
C = {[idx_DORSUM'],[idx_Gyro']};
[uvals, ~, bin] = unique(cell2mat(C));
origarray = repelem(1:numel(C), cellfun(@numel, C));
dist = table(uvals', accumarray(bin, 1), accumarray(bin, origarray, [], @(x) {x}), 'VariableNames', {'index', 'repetition', 'arrayindices'});
dist = sortrows(dist, 'repetition');
idx=find(dist.repetition==2);
idx_DORSUM_Gyro=dist.index(idx);
clear('bin', 'C', 'dist', 'idx', 'origarray','uvals');

% Indices for dorsum accelerometer
C = {[idx_DORSUM'],[idx_Acc']};
[uvals, ~, bin] = unique(cell2mat(C));
origarray = repelem(1:numel(C), cellfun(@numel, C));
dist = table(uvals', accumarray(bin, 1), accumarray(bin, origarray, [], @(x) {x}), 'VariableNames', {'index', 'repetition', 'arrayindices'});
dist = sortrows(dist, 'repetition');
idx=find(dist.repetition==2);
idx_DORSUM_Acc=dist.index(idx);
clear('bin', 'C', 'dist', 'idx', 'origarray','uvals');

% Indices for dorsum magnetometer
C = {[idx_DORSUM'],[idx_Mag']};
[uvals, ~, bin] = unique(cell2mat(C));
origarray = repelem(1:numel(C), cellfun(@numel, C));
dist = table(uvals', accumarray(bin, 1), accumarray(bin, origarray, [], @(x) {x}), 'VariableNames', {'index', 'repetition', 'arrayindices'});
dist = sortrows(dist, 'repetition');
idx=find(dist.repetition==2);
idx_DORSUM_Mag=dist.index(idx);
clear('bin', 'C', 'dist', 'idx', 'origarray','uvals');

%% Find indices for Muscle EMG

% Flexor C Radialis index
C = {[idx_FLXCRAD'],[idx_EMG']};
[uvals, ~, bin] = unique(cell2mat(C));
origarray = repelem(1:numel(C), cellfun(@numel, C));
dist = table(uvals', accumarray(bin, 1), accumarray(bin, origarray, [], @(x) {x}), 'VariableNames', {'index', 'repetition', 'arrayindices'});
dist = sortrows(dist, 'repetition');
idx=find(dist.repetition==2);
idx_FLXCRAD_EMG=dist.index(idx);
clear('bin', 'C', 'dist', 'idx', 'origarray','uvals');

% Extensor C Ulnaris index
C = {[idx_EXTCULN'],[idx_EMG']};
[uvals, ~, bin] = unique(cell2mat(C));
origarray = repelem(1:numel(C), cellfun(@numel, C));
dist = table(uvals', accumarray(bin, 1), accumarray(bin, origarray, [], @(x) {x}), 'VariableNames', {'index', 'repetition', 'arrayindices'});
dist = sortrows(dist, 'repetition');
idx=find(dist.repetition==2);
idx_EXTCULN_EMG=dist.index(idx);
clear('bin', 'C', 'dist', 'idx', 'origarray','uvals');

% Extensor C Radialis index
C = {[idx_EXTCRAD'],[idx_EMG']};
[uvals, ~, bin] = unique(cell2mat(C));
origarray = repelem(1:numel(C), cellfun(@numel, C));
dist = table(uvals', accumarray(bin, 1), accumarray(bin, origarray, [], @(x) {x}), 'VariableNames', {'index', 'repetition', 'arrayindices'});
dist = sortrows(dist, 'repetition');
idx=find(dist.repetition==2);
idx_EXTCRAD_EMG=dist.index(idx);
clear('bin', 'C', 'dist', 'idx', 'origarray','uvals');

%% Find indices for Nerve EMG

% Nerve posterior in EMG index
C = {[idx_NRVPOSTINT'],[idx_EMG']};
[uvals, ~, bin] = unique(cell2mat(C));
origarray = repelem(1:numel(C), cellfun(@numel, C));
dist = table(uvals', accumarray(bin, 1), accumarray(bin, origarray, [], @(x) {x}), 'VariableNames', {'index', 'repetition', 'arrayindices'});
dist = sortrows(dist, 'repetition');
idx=find(dist.repetition==2);
idx_NRVPOSTINT_EMG=dist.index(idx);
clear('bin', 'C', 'dist', 'idx', 'origarray','uvals');

% Nerve posterior out EMG index
C = {[idx_NRVPOSTOUT'],[idx_EMG']};
[uvals, ~, bin] = unique(cell2mat(C));
origarray = repelem(1:numel(C), cellfun(@numel, C));
dist = table(uvals', accumarray(bin, 1), accumarray(bin, origarray, [], @(x) {x}), 'VariableNames', {'index', 'repetition', 'arrayindices'});
dist = sortrows(dist, 'repetition');
idx=find(dist.repetition==2);
idx_NRVPOSTOUT_EMG=dist.index(idx);
clear('bin', 'C', 'dist', 'idx', 'origarray','uvals');

% Nerve anterior in EMG index
C = {[idx_NRVANTINT'],[idx_EMG']};
[uvals, ~, bin] = unique(cell2mat(C));
origarray = repelem(1:numel(C), cellfun(@numel, C));
dist = table(uvals', accumarray(bin, 1), accumarray(bin, origarray, [], @(x) {x}), 'VariableNames', {'index', 'repetition', 'arrayindices'});
dist = sortrows(dist, 'repetition');
idx=find(dist.repetition==2);
idx_NRVANTINT_EMG=dist.index(idx);
clear('bin', 'C', 'dist', 'idx', 'origarray','uvals');

% Nerve anterior out EMG index
C = {[idx_NRVANTOUT'],[idx_EMG']};
[uvals, ~, bin] = unique(cell2mat(C));
origarray = repelem(1:numel(C), cellfun(@numel, C));
dist = table(uvals', accumarray(bin, 1), accumarray(bin, origarray, [], @(x) {x}), 'VariableNames', {'index', 'repetition', 'arrayindices'});
dist = sortrows(dist, 'repetition');
idx=find(dist.repetition==2);
idx_NRVANTOUT_EMG=dist.index(idx);
clear('bin', 'C', 'dist', 'idx', 'origarray','uvals');

%% Remove old indices

clear('idx_Acc','idx_DORSUM','idx_EMG','idx_EXTCRAD','idx_EXTCULN','idx_FLXCRAD','idx_Gyro','idx_Mag',...
    'idx_NRVANTINT','idx_NRVANTOUT','idx_NRVPOSTINT','idx_NRVPOSTOUT');
%% Calculate RMS muscle activity per muscle

% Flexor Carpi Radialis
for idx_var=1:length(zworksp_idx);
    name_var=(zworksp_idx(idx_var).name);
    curr_var=eval(name_var);
    [rms_FLXCRAD ~]=envelope(abs(curr_var(idx_FLXCRAD_EMG,:)-mean(curr_var(idx_FLXCRAD_EMG,:))),400,'rms'); %RMS
    curr_var(idx_FLXCRAD_EMG,:)=rms_FLXCRAD-min(rms_FLXCRAD); %RMS-minimum RMS
    assignin('base',(zworksp_idx(idx_var).name),curr_var);
end
clear('idx_var','curr_var','name_var','rms_FLXCRAD');
% Extensor Carpi Ulnaris
for idx_var=1:length(zworksp_idx);
    name_var=(zworksp_idx(idx_var).name);
    curr_var=eval(name_var);
    [rms_EXTCULN ~]=envelope(abs(curr_var(idx_EXTCULN_EMG,:)-mean(curr_var(idx_EXTCULN_EMG,:))),400,'rms');
    curr_var(idx_EXTCULN_EMG,:)=rms_EXTCULN-min(rms_EXTCULN);
    assignin('base',(zworksp_idx(idx_var).name),curr_var);
end
clear('idx_var','curr_var','name_var','rms_EXTCULN');
% Extensor Carpi Radialis
for idx_var=1:length(zworksp_idx);
    name_var=(zworksp_idx(idx_var).name);
    curr_var=eval(name_var);
    [rms_EXTCRAD ~]=envelope(abs(curr_var(idx_EXTCRAD_EMG,:)-mean(curr_var(idx_EXTCRAD_EMG,:))),400,'rms');
    curr_var(idx_EXTCRAD_EMG,:)=rms_EXTCRAD-min(rms_EXTCRAD);
    assignin('base',(zworksp_idx(idx_var).name),curr_var);
end
clear('idx_var','curr_var','name_var','rms_EXTCRAD');

%% Horizontally concatonate recordings

%all_RadialDeviation=horzcat(RadialDeviation, RadialDeviation_1kg, RadialDeviation_2kg);
%all_WristFlexion=horzcat(WristFlexion, WristFlexion_1kg, WristFlexion_2kg);
all_data=horzcat(RadialDeviation, RadialDeviation_1kg, RadialDeviation_2kg,WristFlexion, WristFlexion_1kg, WristFlexion_2kg);

%% Save individual nerve sensor data for same motions across all recordings

name_var='all_data';
curr_var=eval(name_var);
%
savefile=strcat(name_var,'_','NRVPOSTOUT_EMG.mat');
x=curr_var(idx_NRVPOSTOUT_EMG,:);
save(savefile,'x');
%
savefile=strcat(name_var,'_','NRVPOSTINT_EMG.mat');
x=curr_var(idx_NRVPOSTINT_EMG,:);
save(savefile,'x');
%
savefile=strcat(name_var,'_','NRVANTINT_EMG.mat');
x=curr_var(idx_NRVANTINT_EMG,:);
save(savefile,'x');
%
savefile=strcat(name_var,'_','NRVANTOUT_EMG.mat');
x=curr_var(idx_NRVANTOUT_EMG,:);
save(savefile,'x');
%
all_data_DORSUM_Acc=curr_var(idx_DORSUM_Acc,:);
savefile=strcat(name_var,'_','DORSUM_Acc.mat');
save(savefile,'all_data_DORSUM_Acc');
%
all_data_DORSUM_Gyro=curr_var(idx_DORSUM_Gyro,:);
savefile=strcat(name_var,'_','DORSUM_Gyro.mat');
save(savefile,'all_data_DORSUM_Gyro');
%
all_data_DORSUM_Mag=curr_var(idx_DORSUM_Mag,:);
savefile=strcat(name_var,'_','DORSUM_Mag.mat');
save(savefile,'all_data_DORSUM_Mag');

clear('savefile', 'x', 'name_var', 'curr_var');

%% Save individual nerve EMG data per motion

% for idx_var=1:length(zworksp_idx);
%     name_var=(zworksp_idx(idx_var).name);
%     curr_var=eval(name_var);
%     savefile=strcat(name_var,'_','NRVPOSTOUT_EMG.mat');
%     x=curr_var(idx_NRVPOSTOUT_EMG,:);
%     save(savefile,'x');
%     savefile=strcat(name_var,'_','NRVPOSTINT_EMG.mat');
%     x=curr_var(idx_NRVPOSTINT_EMG,:);
%     save(savefile,'x');
%     savefile=strcat(name_var,'_','NRVANTINT_EMG.mat');
%     x=curr_var(idx_NRVANTINT_EMG,:);
%     save(savefile,'x');
%     savefile=strcat(name_var,'_','NRVANTOUT_EMG.mat');
%     x=curr_var(idx_NRVANTOUT_EMG,:);
%     save(savefile,'x');
% end

%% Run EMG decomposition (LEMG_Analyzer.m)
% Perform this step manually by selecting lines and clicking function key 9

% Enter directory below of editted Sedghamiz(2015) spike sorting toolbox
addpath('C:\Users\twilm\OneDrive\Documents\GitHub\SMC10_ThomasWilmot\SMC10ThomasWilmot\EMG_Analyzer_TWedit\LEMG_Analyzer_V8');

% Open spike sorting toolbox
run LEMG_Analyzer.m

% >open desired EMG or EEG data>enter sample frequency>click start clustering
% 0.05 0.3
display('0.05 0.3 with NRVANTINT works best. Press space to continue.');
pause;
%% Convert Struct of templates into cell

Templates=struct2cell(Templates);

%% Remove MUAP with label '-1'
idx_error=find(MUAP_label==(-1));
MUAP_label(idx_error)=[];
MUAP_idx(idx_error)=[];
clear('idx_error')
%% Save Output from spike sorting to save time later

save('TMW_SMC10_wlemgout.mat')

%% Load output from LEMG_Anaylzer.m
clc;
clear all;
load('TMW_SMC10_wlemgout.mat');
load('all_data_DORSUM_Gyro');
load('all_data_DORSUM_Mag');
load('all_data_DORSUM_Acc');


%% Scores Matrix

% scrs_all(1,:)= Recording Index
% 1= Radial Deviation|none
% 2= Radial Deviation|1kg
% 3= Radial Deviation|2kg
% 4= Wrist Flexion|none
% 5= Wrist Flexion|1kg
% 6= Wrist Flexion|2kg

% scrs_all(2,:)= Motion
% 1= Motion
% 0= Rest

% scrs_all(3,:)= Spike Label during motion

% scrs_all(4,:)= Sample Index

% scrs_all(5,:)= Inter-Spike Interval per spike label
% =number of samples till spike of same label
%% Create labelled time windows from EMG and IMU data (1 of 2)

% Create empty matrix for labelling over time
scrs_all=zeros([5 length(all_data_DORSUM_Acc)]);

wndw=round(30*fs_all)-1; %Recording window

for i=0:5;
    scrs_all((1),(i*wndw)+1:(i*wndw)+wndw)=(i+1)*ones([1,wndw]);
end
clear('i');

%% Sample Index

scrs_all(4,:)=1:(length(scrs_all));

%% Normalise and scale IMU data according to own results

%Acc
%Subtract mean per motion per axis
all_data_DORSUM_Acc(:,1:(3*wndw))=((all_data_DORSUM_Acc(:,...
    1:(3*wndw))')-mean(all_data_DORSUM_Acc(:,1:(3*wndw))'))'; %Rad Dev
all_data_DORSUM_Acc(:,(3*wndw):(6*wndw))=((all_data_DORSUM_Acc(:,...
    (3*wndw):(6*wndw))')-mean(all_data_DORSUM_Acc(:,(3*wndw):(6*wndw))'))'; %Wrst Flx
% Remove recording start and end edges
all_data_DORSUM_Acc(:,(1:1111))=0;
all_data_DORSUM_Acc(:,(end-1111:end))=0;
for i=1:5
    all_data_DORSUM_Acc(:,[(i*wndw)-1111:(i*wndw)+1111])=0;
end
clear('i');
%divide by max per motion
all_data_DORSUM_Acc(:,1:(3*wndw))=all_data_DORSUM_Acc(:,...
    1:(3*wndw))/(max(max(abs(all_data_DORSUM_Acc(:,1:(3*wndw)))))); %Rad Dev
all_data_DORSUM_Acc(:,(3*wndw):(6*wndw))=all_data_DORSUM_Acc(:,...
    (3*wndw):(6*wndw))/(max(max(abs(all_data_DORSUM_Acc(:,(3*wndw):(6*wndw)))))); %Wrst Flx

%Gyro
%Subtract mean per motion per axis
all_data_DORSUM_Gyro(:,1:(3*wndw))=((all_data_DORSUM_Gyro(:,...
    1:(3*wndw))')-mean(all_data_DORSUM_Gyro(:,1:(3*wndw))'))'; %Rad Dev
all_data_DORSUM_Gyro(:,(3*wndw):(6*wndw))=((all_data_DORSUM_Gyro(:,...
    (3*wndw):(6*wndw))')-mean(all_data_DORSUM_Gyro(:,(3*wndw):(6*wndw))'))'; %Wrst Flx
% Remove recording start and end edges
all_data_DORSUM_Gyro(:,(1:1111))=0;
all_data_DORSUM_Gyro(:,(end-1111:end))=0;
for i=1:5
    all_data_DORSUM_Gyro(:,[(i*wndw)-1111:(i*wndw)+1111])=0;
end
clear('i');
%divide by max per motion
all_data_DORSUM_Gyro(:,1:(3*wndw))=all_data_DORSUM_Gyro(:,...
    1:(3*wndw))/(max(max(abs(all_data_DORSUM_Gyro(:,1:(3*wndw)))))); %Rad Dev
all_data_DORSUM_Gyro(:,(3*wndw):(6*wndw))=all_data_DORSUM_Gyro(:,...
    (3*wndw):(6*wndw))/(max(max(abs(all_data_DORSUM_Gyro(:,(3*wndw):(6*wndw)))))); %Wrst Flx

%Mag
%Subtract mean per motion per axis
all_data_DORSUM_Mag(:,1:(3*wndw))=((all_data_DORSUM_Mag(:,...
    1:(3*wndw))')-mean(all_data_DORSUM_Mag(:,1:(3*wndw))'))'; %Rad Dev
all_data_DORSUM_Mag(:,(3*wndw):(6*wndw))=((all_data_DORSUM_Mag(:,...
    (3*wndw):(6*wndw))')-mean(all_data_DORSUM_Mag(:,(3*wndw):(6*wndw))'))'; %Wrst Flx
% Remove recording start and end edges
all_data_DORSUM_Mag(:,(1:1111))=0;
all_data_DORSUM_Mag(:,(end-1111:end))=0;
for i=1:5
    all_data_DORSUM_Mag(:,[(i*wndw)-1111:(i*wndw)+1111])=0;
end
clear('i');
%divide by max per motion
all_data_DORSUM_Mag(:,1:(3*wndw))=all_data_DORSUM_Mag(:,...
    1:(3*wndw))/(max(max(abs(all_data_DORSUM_Mag(:,1:(3*wndw)))))); %Rad Dev
all_data_DORSUM_Mag(:,(3*wndw):(6*wndw))=all_data_DORSUM_Mag(:,...
    (3*wndw):(6*wndw))/(max(max(abs(all_data_DORSUM_Mag(:,(3*wndw):(6*wndw)))))); %Wrst Flx

%% Create labelled time windows from EMG and IMU data (2 of 2)

% Start of contraction (EMG)

% End of contraction (EMG)

% Rest position (minimum muscle activity)

% Held pose position (prolonged muscle activity)

% Onsets and offsets of motion (IMU)
[zpw,zINITCROSS,zFINALCROSS]=pulsewidth(((all_data_DORSUM_Acc(1,...
    :)*(-1)))/max((all_data_DORSUM_Acc(1,:)*(-1))),'StateLevels',[-0.2 0.2]);
idx_error=find(zpw<2222);
zpw(:,idx_error)=[];
zINITCROSS(:,idx_error)=[];
zFINALCROSS(:,idx_error)=[];
zINITCROSS=round(zINITCROSS); %Onset
zFINALCROSS=round(zFINALCROSS); %Offset
clear('idx_error','zpw');
x=[1:length(all_data)];
for i=1:length(zINITCROSS);
    idx_mo{:,i}=find(x > (zINITCROSS(i)) & x < (zFINALCROSS(i)));
end
clear('x');
scrs_all(2,[idx_mo{:,1:length(zINITCROSS)}])=1;

% Percentage muscle strength used (none,1kg,2kg)


%% Classify force applied (low, medium 1kg, higher 2kg)
%% Add spike labels to scores matrix
x=MUAP_idx;
for i=1:length(zINITCROSS);
    idx_mospke{:,i}=find(x>(zINITCROSS(i)) & x<(zFINALCROSS(i))); %index of spike during motion per neuron
end
clear('x');
x=[idx_mospke{1:end}];
for i=1:length(x);
    scrs_all(3,MUAP_idx(x(i)))=MUAP_label(x(i)); %label of spike during motion at idx
end

% Remove samples with no spikes
% idx_clear=find(scrs_all(3,:)==0);
% scrs_all(:,idx_clear)=[];
clear('idx_clear');

%% Calculate Interspike intervals per neuron

for idx=1:max(unique(scrs_all(3,:)));
    idx_neuron{idx}=find(scrs_all(3,:)==idx);
end
clear('idx')

% interspikes intervals for following 5 spikes
for idx=1:max(unique(scrs_all(3,:))); %neuron
    for idx2=1:length(idx_neuron{idx})-5; %spike per neuron
        for idx3=0:4; %range of spikes to check intervals for
            scrs_all(5,(idx_neuron{idx}(idx2))+idx3)=scrs_all(4,(idx_neuron{idx}(idx2+1+idx3)))-scrs_all(4,(idx_neuron{idx}(idx2)));
        end
    end
end

% set spike intervals between motions to 0
idx_across=find(scrs_all(5,:)>2200);
scrs_all(5,idx_across)=0;

%% PCA of scores matrix

% % Motions 1-3
% [coeff1,score1,latent1]=pca(scrs_all(:,1:1128)');%1:(3*wndw)
% figure()
% biplot(coeff1(:,1:2),'scores',score1(:,1:2),'varlabels',{'v1','v2','v3'});
% clear('coeff', 'score', 'latent');
%
% % Motions 4-6
% [coeff2,score2,latent2]=pca(scrs_all(:,1129:end)');%(3*wndw)+1:end
% %figure()
% hold on;
% biplot(coeff2(:,1:2),'scores',score2(:,1:2),'varlabels',{'w1','w2','w3'});
% clear('coeff', 'score', 'latent');

%% Separate scores per motion


%test motion
scrs_tst1=scrs_all(:,(3500):(8000));
scrs_tst2=scrs_all(:,(7.1*10^4):(7.7*10^4));
scrs_tst3=scrs_all(:,(1.85*10^5):(1.92*10^5));
scrs_tst4=scrs_all(:,(2.37*10^5):(2.41*10^5));
scrs_tst5=scrs_all(:,(3*10^5):(3.06*10^5));
scrs_tst6=scrs_all(:,(3.72*10^5):(3.77*10^5));

 %Remove test data from original 
% scrs_all(:,(3.72*10^5):(3.77*10^5))=0;
% scrs_all(:,(3*10^5):(3.06*10^5))=0;
% scrs_all(:,(2.37*10^5):(2.41*10^5))=0;
% scrs_all(:,(1.85*10^5):(1.92*10^5))=0;
% scrs_all(:,(7.1*10^4):(7.7*10^4))=0;
% scrs_all(:,(3500):(8000))=0;

%Recordings without test data segments
scrs_1=scrs_all(:,1:wndw);
scrs_2=scrs_all(:,(wndw+1):(2*wndw));
scrs_3=scrs_all(:,(2*wndw+1):(3*wndw));
scrs_4=scrs_all(:,(3*wndw+1):(4*wndw));
scrs_5=scrs_all(:,(4*wndw+1):(5*wndw));
scrs_6=scrs_all(:,(5*wndw+1):(6*wndw));

%% Remove 0 spike labels

idx=find(scrs_1(3,:)==0);
scrs_1(:,idx)=[];
idx=find(scrs_2(3,:)==0);
scrs_2(:,idx)=[];
idx=find(scrs_3(3,:)==0);
scrs_3(:,idx)=[];
idx=find(scrs_4(3,:)==0);
scrs_4(:,idx)=[];
idx=find(scrs_5(3,:)==0);
scrs_5(:,idx)=[];
idx=find(scrs_6(3,:)==0);
scrs_6(:,idx)=[];

%test motions
idx=find(scrs_tst1(3,:)==0);
scrs_tst1(:,idx)=[];
idx=find(scrs_tst2(3,:)==0);
scrs_tst2(:,idx)=[];
idx=find(scrs_tst3(3,:)==0);
scrs_tst3(:,idx)=[];
idx=find(scrs_tst4(3,:)==0);
scrs_tst4(:,idx)=[];
idx=find(scrs_tst5(3,:)==0);
scrs_tst5(:,idx)=[];
idx=find(scrs_tst6(3,:)==0);
scrs_tst6(:,idx)=[];

%% Remove (recording) and (binary motion)
scrs_1(1:2,:)=[];
scrs_2(1:2,:)=[];
scrs_3(1:2,:)=[];
scrs_4(1:2,:)=[];
scrs_5(1:2,:)=[];
scrs_6(1:2,:)=[];
%test motion
scrs_tst1(1:2,:)=[];
scrs_tst2(1:2,:)=[];
scrs_tst3(1:2,:)=[];
scrs_tst4(1:2,:)=[];
scrs_tst5(1:2,:)=[];
scrs_tst6(1:2,:)=[];

%% Remove interspike intervals = 0

idx=find(scrs_1(3,:)==0);
scrs_1(:,idx)=[];
idx=find(scrs_2(3,:)==0);
scrs_2(:,idx)=[];
idx=find(scrs_3(3,:)==0);
scrs_3(:,idx)=[];
idx=find(scrs_4(3,:)==0);
scrs_4(:,idx)=[];
idx=find(scrs_5(3,:)==0);
scrs_5(:,idx)=[];
idx=find(scrs_6(3,:)==0);
scrs_6(:,idx)=[];

%test motions
idx=find(scrs_tst1(3,:)==0);
scrs_tst1(:,idx)=[];
idx=find(scrs_tst2(3,:)==0);
scrs_tst2(:,idx)=[];
idx=find(scrs_tst3(3,:)==0);
scrs_tst3(:,idx)=[];
idx=find(scrs_tst4(3,:)==0);
scrs_tst4(:,idx)=[];
idx=find(scrs_tst5(3,:)==0);
scrs_tst5(:,idx)=[];
idx=find(scrs_tst6(3,:)==0);
scrs_tst6(:,idx)=[];

%% Create histograms per neuron per recording
clear('idx2');

for idx=1:max(unique(scrs_2(1,:)));
    %scrs1
    idx2{idx}=find(scrs_1(1,:)==idx);
    [No_ISI1{idx} ~]=histcounts(scrs_1(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot(smooth(No_ISI1{idx},22));hold on; %11 samples=10ms window
    %scrs2
    idx2{idx}=find(scrs_2(1,:)==idx);
    [No_ISI2{idx} ~]=histcounts(scrs_2(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot(smooth(No_ISI2{idx},11));hold on;
    %scrs3
    idx2{idx}=find(scrs_3(1,:)==idx);
    [No_ISI3{idx} ~]=histcounts(scrs_3(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot(smooth(No_ISI3{idx},11),'r--');hold on;
    %scrs4
    idx2{idx}=find(scrs_4(1,:)==idx);
    [No_ISI4{idx} ~]=histcounts(scrs_4(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot(smooth(No_ISI4{idx},11));hold on; %11 samples=5ms window
    %scrs5
    idx2{idx}=find(scrs_5(1,:)==idx);
    [No_ISI5{idx} ~]=histcounts(scrs_5(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot(smooth(No_ISI5{idx},11));hold on;
    %scrs6
    idx2{idx}=find(scrs_6(1,:)==idx);
    [No_ISI6{idx} ~]=histcounts(scrs_6(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot(smooth(No_ISI6{idx},11));hold on;  
    %test motions
    %scrs_tst1
    idx2{idx}=find(scrs_tst1(1,:)==idx);
    [No_ISI_tst1{idx} ~]=histcounts(scrs_tst1(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot(smooth(No_ISI_tst1{idx},22));hold on; %11 samples=10ms window
    %scrs_tst2
    idx2{idx}=find(scrs_tst2(1,:)==idx);
    [No_ISI_tst2{idx} ~]=histcounts(scrs_tst2(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot(smooth(No_ISI_tst2{idx},11));hold on;
    %scrs_tst3
    idx2{idx}=find(scrs_tst3(1,:)==idx);
    [No_ISI_tst3{idx} ~]=histcounts(scrs_tst3(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot(smooth(No_ISI_tst3{idx},11),'r--');hold on;
    %scrs_tst4
    idx2{idx}=find(scrs_tst4(1,:)==idx);
    [No_ISI_tst4{idx} ~]=histcounts(scrs_tst4(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot(smooth(No_ISI_tst4{idx},11));hold on; %11 samples=5ms window
    %scrs_tst5
    idx2{idx}=find(scrs_tst5(1,:)==idx);
    [No_ISI_tst5{idx} ~]=histcounts(scrs_tst5(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot(smooth(No_ISI_tst5{idx},11));hold on;
    %scrs_tst6
    idx2{idx}=find(scrs_tst6(1,:)==idx);
    [No_ISI_tst6{idx} ~]=histcounts(scrs_tst6(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot(smooth(No_ISI_tst6{idx},11));hold on;
end

for idx=1:max(unique(scrs_2(1,:)));
    %scrs3
    idx2{idx}=find(scrs_3(1,:)==idx);
   [No_ISI3{idx} ~]=histcounts(scrs_3(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot((smooth(No_ISI3{idx},11))/max(smooth(No_ISI3{idx},11)),'r--');hold on;
    %scrs_tst
    idx2{idx}=find(scrs_tst3(1,:)==idx);
    [No_ISI_tst3{idx} ~]=histcounts(scrs_tst3(3,idx2{idx}),'BinWidth',1,'BinLimits',[5,1000]);
    figure(idx);
    plot((smooth(No_ISI_tst3{idx},11))/max(smooth(No_ISI_tst3{idx},11)),'k');hold on;
end 
%add all intervals
%smooth/parsing acrtoss probability
%compare autocorrelation
%% Denoise Rate templates

%Check correlation of rate templates to other known rate templates and
%weight by number of spikes per neuron. This gives the ratio of similarity
%between templates
for idx=1:length(idx_neuron);
    %1
    ISImat_all=[No_ISI1(idx) No_ISI2(idx) No_ISI3(idx) No_ISI4(idx) No_ISI5(idx) No_ISI6(idx)]';
    ISImat_all=cell2mat(ISImat_all);

    % Correlation between templates weighted by number of spikes per neuron
    % per motion
    CrTemp2Temp=corr(ISImat_all')*length(find(cell2mat(No_ISI1(idx))>0));
    CrTemp_1(idx,:)=CrTemp2Temp(1,:);
    
    CrTemp2Temp=corr(ISImat_all')*length(find(cell2mat(No_ISI2(idx))>0));
    CrTemp_2(idx,:)=CrTemp2Temp(2,:);
    
    CrTemp2Temp=corr(ISImat_all')*length(find(cell2mat(No_ISI3(idx))>0));
    CrTemp_3(idx,:)=CrTemp2Temp(3,:);
    
    CrTemp2Temp=corr(ISImat_all')*length(find(cell2mat(No_ISI4(idx))>0));
    CrTemp_4(idx,:)=CrTemp2Temp(4,:);
    
    CrTemp2Temp=corr(ISImat_all')*length(find(cell2mat(No_ISI5(idx))>0));
    CrTemp_5(idx,:)=CrTemp2Temp(5,:);
    
    CrTemp2Temp=corr(ISImat_all')*length(find(cell2mat(No_ISI6(idx))>0));
    CrTemp_6(idx,:)=CrTemp2Temp(6,:);
end
clear('CrTemp2Temp')

%Normalise similarity matrix
CrTemp_1=CrTemp_1./max(CrTemp_1')';
CrTemp_2=CrTemp_2./max(CrTemp_2')';
CrTemp_3=CrTemp_3./max(CrTemp_3')';
CrTemp_4=CrTemp_4./max(CrTemp_4')';
CrTemp_5=CrTemp_5./max(CrTemp_5')';
CrTemp_6=CrTemp_6./max(CrTemp_6')';

%Remove negative
idx=find(CrTemp_1<0);
CrTemp_1(idx)=0;
idx=find(CrTemp_2<0);
CrTemp_2(idx)=0;
idx=find(CrTemp_3<0);
CrTemp_3(idx)=0;
idx=find(CrTemp_4<0);
CrTemp_4(idx)=0;
idx=find(CrTemp_5<0);
CrTemp_5(idx)=0;
idx=find(CrTemp_6<0);
CrTemp_6(idx)=0;


%Remove ratio of similar templates
for idx=1:length(idx_neuron)
    %1
    No_ISI1n(idx,:)=((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))-...
        nansum([((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_1(idx,2)));...
        ((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_1(idx,3)));...
        ((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_1(idx,4)));...
        ((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_1(idx,5)));...
        ((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_1(idx,6)))]);
    %2
     No_ISI2n(idx,:)=((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))-...
        nansum([((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_2(idx,1)));...
        ((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_2(idx,3)));...
        ((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_2(idx,4)));...
        ((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_2(idx,5)));...
        ((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_2(idx,6)))]);
    %3
     No_ISI3n(idx,:)=((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))-...
        nansum([((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_3(idx,2)));...
        ((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_3(idx,1)));...
        ((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_3(idx,4)));...
        ((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_3(idx,5)));...
        ((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_3(idx,6)))]);
    %4
    No_ISI4n(idx,:)=((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))-...
        nansum([((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_4(idx,2)));...
        ((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_4(idx,3)));...
        ((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_4(idx,1)));...
        ((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_4(idx,5)));...
        ((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_4(idx,6)))]);
    %5
 No_ISI5n(idx,:)=((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))-...
        nansum([((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_5(idx,2)));...
        ((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_5(idx,3)));...
        ((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_5(idx,4)));...
        ((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_5(idx,1)));...
        ((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_5(idx,6)))]);
%6
     No_ISI6n(idx,:)=((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))-...
        nansum([((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_6(idx,2)));...
        ((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_6(idx,3)));...
        ((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_6(idx,4)));...
        ((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_6(idx,5)));...
        ((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_6(idx,1)))]);
end
%Remove negative values
idx=find(No_ISI1n<0);
No_ISI1n(idx)=0;
idx=find(No_ISI2n<0);
No_ISI2n(idx)=0;
idx=find(No_ISI3n<0);
No_ISI3n(idx)=0;
idx=find(No_ISI4n<0);
No_ISI4n(idx)=0;
idx=find(No_ISI5n<0);
No_ISI5n(idx)=0;
idx=find(No_ISI6n<0);
No_ISI6n(idx)=0;

%clear('No_ISI1', 'No_ISI2', 'No_ISI3', 'No_ISI4', 'No_ISI5', 'No_ISI6');
%% Correlation with test motion and motion templates
for idx=1:length(idx_neuron);
%tst1
ISImat_tst1=[No_ISI1n(idx,:) No_ISI2n(idx,:) No_ISI3n(idx,:) No_ISI4n(idx,:) No_ISI5n(idx,:) No_ISI6n(idx,:) No_ISI_tst1(idx)]';
ISImat_tst1=cell2mat(ISImat_tst1);
%tst2
ISImat_tst2=[No_ISI1n(idx,:) No_ISI2n(idx,:) No_ISI3n(idx,:) No_ISI4n(idx,:) No_ISI5n(idx,:) No_ISI6n(idx,:) No_ISI_tst2(idx)]';
ISImat_tst2=cell2mat(ISImat_tst2);
%tst3
ISImat_tst3=[No_ISI1n(idx,:) No_ISI2n(idx,:) No_ISI3n(idx,:) No_ISI4n(idx,:) No_ISI5n(idx,:) No_ISI6n(idx,:) No_ISI_tst3(idx)]';
ISImat_tst3=cell2mat(ISImat_tst3);
%tst4
ISImat_tst4=[No_ISI1n(idx,:) No_ISI2n(idx,:) No_ISI3n(idx,:) No_ISI4n(idx,:) No_ISI5n(idx,:) No_ISI6n(idx,:) No_ISI_tst4(idx)]';
ISImat_tst4=cell2mat(ISImat_tst4);
%tst5
ISImat_tst5=[No_ISI1n(idx,:) No_ISI2n(idx,:) No_ISI3n(idx,:) No_ISI4n(idx,:) No_ISI5n(idx,:) No_ISI6n(idx,:) No_ISI_tst5(idx)]';
ISImat_tst5=cell2mat(ISImat_tst5);
%tst6
ISImat_tst6=[No_ISI1n(idx,:) No_ISI2n(idx,:) No_ISI3n(idx,:) No_ISI4n(idx,:) No_ISI5n(idx,:) No_ISI6n(idx,:) No_ISI_tst6(idx)]';
ISImat_tst6=cell2mat(ISImat_tst6);

MLM=corr(ISImat_tst1')*length(find(ISImat_tst1>0));%Most Likely Motion
MLMtst_1(idx,:)=MLM(7,:);

MLM=corr(ISImat_tst2')*length(find(ISImat_tst2>0));%Most Likely Motion
MLMtst_2(idx,:)=MLM(7,:);

MLM=corr(ISImat_tst3')*length(find(ISImat_tst3>0));%Most Likely Motion
MLMtst_3(idx,:)=MLM(7,:);

MLM=corr(ISImat_tst4')*length(find(ISImat_tst4>0));%Most Likely Motion
MLMtst_4(idx,:)=MLM(7,:);

MLM=corr(ISImat_tst5')*length(find(ISImat_tst5>0));%Most Likely Motion
MLMtst_5(idx,:)=MLM(7,:);

MLM=corr(ISImat_tst6')*length(find(ISImat_tst6>0));%Most Likely Motion
MLMtst_6(idx,:)=MLM(7,:);
end
clear('MLM')
figure(21);
bar(((nansum(MLMtst_1)))/max((nansum(MLMtst_1)))); hold on;
ylabel('Correlation');
xlabel('Motion Index (7 is test motion)');
title('Motion 1: Radial Deviation 0Kg');
figure(22);
bar(((nansum(MLMtst_2)))/max((nansum(MLMtst_2)))); hold on;
ylabel('Correlation');
xlabel('Motion Index (7 is test motion)');
title('Motion 2: Radial Deviation 1Kg');
figure(23);
bar(((nansum(MLMtst_3)))/max((nansum(MLMtst_3)))); hold on;
ylabel('Correlation');
xlabel('Motion Index (7 is test motion)');
title('Motion 3: Radial Deviation 2Kg');
figure(24);
bar(((nansum(MLMtst_4)))/max((nansum(MLMtst_4)))); hold on;
ylabel('Correlation');
xlabel('Motion Index (7 is test motion)');
title('Motion 4: Wrist Flexion 0Kg');
figure(25);
bar(((nansum(MLMtst_5)))/max((nansum(MLMtst_5)))); hold on;
ylabel('Correlation');
xlabel('Motion Index (7 is test motion)');
title('Motion 5: Wrist Flexion 1Kg');
figure(26);
bar(((nansum(MLMtst_6)))/max((nansum(MLMtst_6)))); hold on;
ylabel('Correlation');
xlabel('Motion Index (7 is test motion)');
title('Motion 6: Wrist Flexion 2Kg');
%% Check
% for idx=1:length(idx_neuron);
% %1
% ISImat_all=[No_ISI1n(idx,:) No_ISI2n(idx,:) No_ISI3n(idx,:) No_ISI4n(idx,:) No_ISI5n(idx,:) No_ISI6n(idx,:) No_ISInew(idx,:)]';
% ISImat_all=cell2mat(ISImat_all);
%
%
% MLM=corr(ISImat_all')*length(find((No_ISInew)>0));%Most Likely Motion
% MLM_new(idx,:)=MLM(7,:);
%
% end
% clear('MLM')
% MLM_new(:,7)=[];
% %Normalise similarity matrix
% MLM_new=MLM_new./max(MLM_new')';
%
% %Remove negative
% idx=find(MLM_new<0);
% MLM_new(idx)=0;

%% Spike rate per neuron per motion weighted by contraction time
% for i=1:length(zINITCROSS);
%  for i2=1:length(idx_lbl);
%      spike_rate(i2,i)=(length(find((scrs_all(3,...
%          zINITCROSS(i):zFINALCROSS(i)))==...
%          i2)))/length(scrs_all(3,zINITCROSS(i):zFINALCROSS(i)));
%  end
%  end


%% Classify primary muscle used

%% Classify directedness/smoothness

%% Score spikes within time windows according to muscle and motion data

%% Compare spike counts from similar motions over increased weight applied during motion

%Repeat below for each motion
% C={[MUAP_label]};
% [uvals,~,bin]=unique(cell2mat(C));
% origarray=repelem(1:numel(C),cellfun(@numel,C));
% dist=table(uvals',accumarray(bin,1),accumarray(bin,origarray,[],@(x) {x}),'VariableNames',{'spikeID_label','repetition','arrayindices'});
% dist=sortrows(dist,'repetition');
%%  Create cells containing onsets, offsets and duration of contractions

% 1=FLXCRAD
% for idx_var=1:length(zworksp_idx);
% name_var=(zworksp_idx(idx_var).name);
% curr_var=eval(name_var);
%         [zpw{idx_var,1},zINITCROSS{idx_var,1},...
%             zFINALCROSS{idx_var,1}]=pulsewidth((curr_var(idx_FLXCRAD_EMG,:))/max(curr_var(idx_FLXCRAD_EMG,:)));
% end

% 2=EXTCRAD
% for idx_var=1:length(zworksp_idx);
% name_var=(zworksp_idx(idx_var).name);
% curr_var=eval(name_var);
%         [zpw{idx_var,2},zINITCROSS{idx_var,2},...
%             zFINALCROSS{idx_var,2}]=pulsewidth((curr_var(idx_EXTCRAD_EMG,:))/max(curr_var(idx_EXTCRAD_EMG,:)));
% end

% 3=EXTCULN
% for idx_var=1:length(zworksp_idx);
% name_var=(zworksp_idx(idx_var).name);
% curr_var=eval(name_var);
%         [zpw{idx_var,3},zINITCROSS{idx_var,3},...
%             zFINALCROSS{idx_var,3}]=pulsewidth((curr_var(idx_EXTCULN_EMG,:))/max(curr_var(idx_EXTCULN_EMG,:)));
% end
% clear( 'curr_var','idx_var','idx_muscle', 'name_var');

%% Total contractile force applied per muscle per motion

% for idx_var=1:length(zworksp_idx);
% name_var=(zworksp_idx(idx_var).name);
% curr_var=eval(name_var);
% % Total force over entire recording per muscle
% totalforce_FLXCRAD(idx_var)=max(cumsum(curr_var(idx_FLXCRAD_EMG,:)));
% totalforce_EXTCRAD(idx_var)=max(cumsum(curr_var(idx_EXTCRAD_EMG,:)));
% totalforce_EXTCULN(idx_var)=max(cumsum(curr_var(idx_EXTCULN_EMG,:)));
% end
% clear('idx_var','curr_var','name_var');
