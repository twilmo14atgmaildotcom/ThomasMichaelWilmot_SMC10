%% Guitar Playing and Motor Interaction
% Thomas Michael Wilmot
% Student of Sound and Music Computing MSc at Aalborg Univeristy
% Last updated: 25/09/18
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
%cd('C:\Users\twilm\OneDrive\Documents\GitHub\SMC10_ThomasWilmot\SMC10ThomasWilmot\Delsysdata210318\TomSimpleMotionDelsys_210318');

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
% 0.05 0.32
display('0.05 0.2 with NRVANTINT works best. Press space to continue.');
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
close all;
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

% Muscle Activty+Motion Onsets/Offsets indices
for idx=1:length(idx_mo);
    OnOff(1,idx)=min(idx_mo{idx});
    OnOff(2,idx)=max(idx_mo{idx});
    OnOff(3,idx)=scrs_all(1,OnOff(1,idx));
end
OnOff(4,:)=OnOff(2,:)-OnOff(1,:);
OnOff(4,:)=OnOff(4,:)/fs_all;
avg_cont=median(OnOff(4,:)); %Median Contraction length
std_cont=std(OnOff(4,:)); %Variance in Contraction Length
idx_nf1=find(OnOff(3,:)==1);
idx_nf2=find(OnOff(3,:)==2);
idx_nf3=find(OnOff(3,:)==3);
idx_nf4=find(OnOff(3,:)==4);
idx_nf5=find(OnOff(3,:)==5);
idx_nf6=find(OnOff(3,:)==6);

for idx_ZTST=1:5;
    %test Motions
    scrs_tst1=scrs_all(:,(OnOff(1,idx_nf1(idx_ZTST))):(OnOff(2,idx_nf1(idx_ZTST))));
    scrs_tst2=scrs_all(:,(OnOff(1,idx_nf2(idx_ZTST))):(OnOff(2,idx_nf2(idx_ZTST))));
    scrs_tst3=scrs_all(:,(OnOff(1,idx_nf3(idx_ZTST))):(OnOff(2,idx_nf3(idx_ZTST))));
    scrs_tst4=scrs_all(:,(OnOff(1,idx_nf4(idx_ZTST))):(OnOff(2,idx_nf4(idx_ZTST))));
    scrs_tst5=scrs_all(:,(OnOff(1,idx_nf5(idx_ZTST))):(OnOff(2,idx_nf5(idx_ZTST))));
    scrs_tst6=scrs_all(:,(OnOff(1,idx_nf6(idx_ZTST))):(OnOff(2,idx_nf6(idx_ZTST))));
    
    scrs_1=[scrs_all(:,1:OnOff(1,idx_nf1(idx_ZTST))),scrs_all(:,OnOff(2,idx_nf1(idx_ZTST)):wndw)];
    scrs_2=[scrs_all(:,wndw+1:OnOff(1,idx_nf2(idx_ZTST))),scrs_all(:,OnOff(2,idx_nf2(idx_ZTST)):2*wndw)];
    scrs_3=[scrs_all(:,2*wndw+1:OnOff(1,idx_nf3(idx_ZTST))),scrs_all(:,OnOff(2,idx_nf3(idx_ZTST)):3*wndw)];
    scrs_4=[scrs_all(:,3*wndw+1:OnOff(1,idx_nf4(idx_ZTST))),scrs_all(:,OnOff(2,idx_nf4(idx_ZTST)):4*wndw)];
    scrs_5=[scrs_all(:,4*wndw+1:OnOff(1,idx_nf5(idx_ZTST))),scrs_all(:,OnOff(2,idx_nf5(idx_ZTST)):5*wndw)];
    scrs_6=[scrs_all(:,5*wndw+1:OnOff(1,idx_nf6(idx_ZTST))),scrs_all(:,OnOff(2,idx_nf6(idx_ZTST)):6*wndw)];
    
    %     training data
    %     scrs_1=scrs_all(:,1:wndw);
    %     scrs_2=scrs_all(:,(wndw+1):(2*wndw));
    %     scrs_3=scrs_all(:,(2*wndw+1):(3*wndw));
    %     scrs_4=scrs_all(:,(3*wndw+1):(4*wndw));
    %     scrs_5=scrs_all(:,(4*wndw+1):(5*wndw));
    %     scrs_6=scrs_all(:,(5*wndw+1):(6*wndw));
    
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
        [No_ISI1{idx} ~]=histcounts((scrs_1(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot(smooth(No_ISI1{idx},5e-03));hold on;
        %scrs2
        idx2{idx}=find(scrs_2(1,:)==idx);
        [No_ISI2{idx} ~]=histcounts((scrs_2(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot(smooth(No_ISI2{idx},5e-03));hold on;
        %scrs3
        idx2{idx}=find(scrs_3(1,:)==idx);
        [No_ISI3{idx} ~]=histcounts((scrs_3(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot(smooth(No_ISI3{idx},5e-03),'r--');hold on;
        %scrs4
        idx2{idx}=find(scrs_4(1,:)==idx);
        [No_ISI4{idx} ~]=histcounts((scrs_4(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot(smooth(No_ISI4{idx},5e-03));hold on; %11 samples=5ms window
        %scrs5
        idx2{idx}=find(scrs_5(1,:)==idx);
        [No_ISI5{idx} ~]=histcounts((scrs_5(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot(smooth(No_ISI5{idx},5e-03));hold on;
        %scrs6
        idx2{idx}=find(scrs_6(1,:)==idx);
        [No_ISI6{idx} ~]=histcounts((scrs_6(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot(smooth(No_ISI6{idx},5e-03));hold on;
        %test motions
        %scrs_tst1
        idx2{idx}=find(scrs_tst1(1,:)==idx);
        [No_ISI_tst1{idx} ~]=histcounts((scrs_tst1(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot(smooth(No_ISI_tst1{idx},22));hold on; %11 samples=10ms window
        %scrs_tst2
        idx2{idx}=find(scrs_tst2(1,:)==idx);
        [No_ISI_tst2{idx} ~]=histcounts((scrs_tst2(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot(smooth(No_ISI_tst2{idx},5e-03));hold on;
        %scrs_tst3
        idx2{idx}=find(scrs_tst3(1,:)==idx);
        [No_ISI_tst3{idx} ~]=histcounts((scrs_tst3(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot(smooth(No_ISI_tst3{idx},5e-03),'r--');hold on;
        %scrs_tst4
        idx2{idx}=find(scrs_tst4(1,:)==idx);
        [No_ISI_tst4{idx} ~]=histcounts((scrs_tst4(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot(smooth(No_ISI_tst4{idx},5e-03));hold on; %11 samples=5ms window
        %scrs_tst5
        idx2{idx}=find(scrs_tst5(1,:)==idx);
        [No_ISI_tst5{idx} ~]=histcounts((scrs_tst5(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot(smooth(No_ISI_tst5{idx},5e-03));hold on;
        %scrs_tst6
        idx2{idx}=find(scrs_tst6(1,:)==idx);
        [No_ISI_tst6{idx} ~]=histcounts((scrs_tst6(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot(smooth(No_ISI_tst6{idx},5e-03));hold on;
    end
    
    for idx=1:max(unique(scrs_2(1,:)));
        %scrs3
        idx2{idx}=find(scrs_3(1,:)==idx);
        [No_ISI3{idx} ~]=histcounts((scrs_3(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot((smooth(No_ISI3{idx},5e-03))/max(smooth(No_ISI3{idx},5e-03)),'r--');hold on;
        %scrs_tst
        idx2{idx}=find(scrs_tst3(1,:)==idx);
        [No_ISI_tst3{idx} ~]=histcounts((scrs_tst3(3,idx2{idx})/fs_all),'BinWidth',1e-03,'BinLimits',[5e-03,140e-03]);
        %figure(idx);
        %plot((smooth(No_ISI_tst3{idx},5e-03))/max(smooth(No_ISI_tst3{idx},5e-03)),'k');hold on;
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
        %         %1
        %         No_ISI1n(idx,:)=(cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx)));
        %         %2
        %         No_ISI2n(idx,:)=(cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx)));
        %
        %         %3
        %         No_ISI3n(idx,:)=(cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx)));
        %         %4
        %         No_ISI4n(idx,:)=(cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx)));
        %
        %         %5
        %         No_ISI5n(idx,:)=(cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx)));
        %         %6
        %         No_ISI6n(idx,:)=(cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx)));
        %         %         %1
        %                 No_ISI1n(idx,:)=((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))-...
        %                     nansum([((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_1(idx,2)));...
        %                     ((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_1(idx,3)));...
        %                     ((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_1(idx,4)));...
        %                     ((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_1(idx,5)));...
        %                     ((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_1(idx,6)))]);
        %         %         %2
        %                 No_ISI2n(idx,:)=((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))-...
        %                     nansum([((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_2(idx,1)));...
        %                     ((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_2(idx,3)));...
        %                     ((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_2(idx,4)));...
        %                     ((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_2(idx,5)));...
        %                     ((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_2(idx,6)))]);
        %         %         %3
        %                 No_ISI3n(idx,:)=((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))-...
        %                     nansum([((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_3(idx,2)));...
        %                     ((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_3(idx,1)));...
        %                     ((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_3(idx,4)));...
        %                     ((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_3(idx,5)));...
        %                     ((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_3(idx,6)))]);
        %         %         %4
        %                 No_ISI4n(idx,:)=((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))-...
        %                     nansum([((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_4(idx,2)));...
        %                     ((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_4(idx,3)));...
        %                     ((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_4(idx,1)));...
        %                     ((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_4(idx,5)));...
        %                     ((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_4(idx,6)))]);
        %         %         %5
        %                 No_ISI5n(idx,:)=((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))-...
        %                     nansum([((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_5(idx,2)));...
        %                     ((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_5(idx,3)));...
        %                     ((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_5(idx,4)));...
        %                     ((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_5(idx,1)));...
        %                     ((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_5(idx,6)))]);
        %         %         %6
        %                 No_ISI6n(idx,:)=((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))-...
        %                     nansum([((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_6(idx,2)));...
        %                     ((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_6(idx,3)));...
        %                     ((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_6(idx,4)));...
        %                     ((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_6(idx,5)));...
        %                     ((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_6(idx,1)))]);
        
        
        
        %1
        No_ISI1n(idx,:)=((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))-...
            nansum([-((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_1(idx,2)));...
            -((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_1(idx,3)));...
            ((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_1(idx,4)));...
            ((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_1(idx,5)));...
            ((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_1(idx,6)))]);
        %2
        No_ISI2n(idx,:)=((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))-...
            nansum([-((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_2(idx,1)));...
            -((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_2(idx,3)));...
            ((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_2(idx,4)));...
            ((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_2(idx,5)));...
            ((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_2(idx,6)))]);
        %3
        No_ISI3n(idx,:)=((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))-...
            nansum([-((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_3(idx,2)));...
            -((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_3(idx,1)));...
            ((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_3(idx,4)));...
            ((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_3(idx,5)));...
            ((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_3(idx,6)))]);
        %4
        No_ISI4n(idx,:)=((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))-...
            nansum([((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_4(idx,2)));...
            ((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_4(idx,3)));...
            ((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_4(idx,1)));...
            -((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_4(idx,5)));...
            -((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_4(idx,6)))]);
        %5
        No_ISI5n(idx,:)=((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))-...
            nansum([((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_5(idx,2)));...
            ((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_5(idx,3)));...
            -((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_5(idx,4)));...
            ((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_5(idx,1)));...
            -((((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))*CrTemp_5(idx,6)))]);
        %6
        No_ISI6n(idx,:)=((cell2mat(No_ISI6(idx)))/max(cell2mat(No_ISI6(idx))))-...
            nansum([((((cell2mat(No_ISI2(idx)))/max(cell2mat(No_ISI2(idx))))*CrTemp_6(idx,2)));...
            ((((cell2mat(No_ISI3(idx)))/max(cell2mat(No_ISI3(idx))))*CrTemp_6(idx,3)));...
            -((((cell2mat(No_ISI4(idx)))/max(cell2mat(No_ISI4(idx))))*CrTemp_6(idx,4)));...
            -((((cell2mat(No_ISI5(idx)))/max(cell2mat(No_ISI5(idx))))*CrTemp_6(idx,5)));...
            ((((cell2mat(No_ISI1(idx)))/max(cell2mat(No_ISI1(idx))))*CrTemp_6(idx,1)))]);
    end
    
    %     end
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
    %normalise test hiscounts
    % for idx=1:length(No_ISI_tst1);
    % No_ISI_tst1{idx}=No_ISI_tst1{idx}/nanmax(No_ISI_tst1{idx}');
    % No_ISI_tst2{idx}=No_ISI_tst2{idx}/nanmax(No_ISI_tst2{idx}');
    % No_ISI_tst3{idx}=No_ISI_tst3{idx}/nanmax(No_ISI_tst3{idx}');
    % No_ISI_tst4{idx}=No_ISI_tst4{idx}/nanmax(No_ISI_tst4{idx}');
    % No_ISI_tst5{idx}=No_ISI_tst5{idx}/nanmax(No_ISI_tst5{idx}');
    % No_ISI_tst6{idx}=No_ISI_tst6{idx}/nanmax(No_ISI_tst6{idx}');
    % end
    clear('No_ISI1', 'No_ISI2', 'No_ISI3', 'No_ISI4', 'No_ISI5', 'No_ISI6');
    
    
    %% Correlation with test motion and motion templates
    %     for idx=1:length(idx_neuron);
    %         %tst1
    %         ISImat_tst1=[No_ISI1n(idx,:)  No_ISI4n(idx,:) No_ISI_tst1(idx)]';
    %         ISImat_tst1=cell2mat(ISImat_tst1);
    %         %tst2
    %         ISImat_tst2=[No_ISI1n(idx,:)  No_ISI4n(idx,:) No_ISI_tst2(idx)]';
    %         ISImat_tst2=cell2mat(ISImat_tst2);
    %         %tst3
    %         ISImat_tst3=[No_ISI1n(idx,:)  No_ISI4n(idx,:) No_ISI_tst3(idx)]';
    %         ISImat_tst3=cell2mat(ISImat_tst3);
    %         %tst4
    %         ISImat_tst4=[No_ISI1n(idx,:)  No_ISI4n(idx,:) No_ISI_tst4(idx)]';
    %         ISImat_tst4=cell2mat(ISImat_tst4);
    %         %tst5
    %         ISImat_tst5=[No_ISI1n(idx,:)  No_ISI4n(idx,:) No_ISI_tst5(idx)]';
    %         ISImat_tst5=cell2mat(ISImat_tst5);
    %         %tst6
    %         ISImat_tst6=[No_ISI1n(idx,:)  No_ISI4n(idx,:) No_ISI_tst6(idx)]';
    %         ISImat_tst6=cell2mat(ISImat_tst6);
    %
    %         MLM=corr(ISImat_tst1')*length(find(ISImat_tst1>0));%Most Likely Motion
    %         MLMtst_1(idx,:)=MLM(7,:);
    %
    %         MLM=corr(ISImat_tst2')*length(find(ISImat_tst2>0));%Most Likely Motion
    %         MLMtst_2(idx,:)=MLM(7,:);
    %
    %         MLM=corr(ISImat_tst3')*length(find(ISImat_tst3>0));%Most Likely Motion
    %         MLMtst_3(idx,:)=MLM(7,:);
    %
    %         MLM=corr(ISImat_tst4')*length(find(ISImat_tst4>0));%Most Likely Motion
    %         MLMtst_4(idx,:)=MLM(7,:);
    %
    %         MLM=corr(ISImat_tst5')*length(find(ISImat_tst5>0));%Most Likely Motion
    %         MLMtst_5(idx,:)=MLM(7,:);
    %
    %         MLM=corr(ISImat_tst6')*length(find(ISImat_tst6>0));%Most Likely Motion
    %         MLMtst_6(idx,:)=MLM(7,:);
    %     end
    %     clear('MLM')
    %     RESULTS1(idx_ZTST,:)=((nansum(MLMtst_1)))/max((nansum(MLMtst_1)));
    %     RESULTS2(idx_ZTST,:)=((nansum(MLMtst_2)))/max((nansum(MLMtst_2)));
    %     RESULTS3(idx_ZTST,:)=((nansum(MLMtst_3)))/max((nansum(MLMtst_3)));
    %     RESULTS4(idx_ZTST,:)=((nansum(MLMtst_4)))/max((nansum(MLMtst_4)));
    %     RESULTS5(idx_ZTST,:)=((nansum(MLMtst_5)))/max((nansum(MLMtst_5)));
    %     RESULTS6(idx_ZTST,:)=((nansum(MLMtst_6)))/max((nansum(MLMtst_6)));
    % end
    %  for idx=1:length(idx_neuron);
    % %         %tst1
    % %         ISImat_tst1=[sum([No_ISI1n(idx,:); No_ISI2n(idx,:); No_ISI3n(idx,:)]) sum([No_ISI4n(idx,:); No_ISI5n(idx,:); No_ISI6n(idx,:)]) No_ISI_tst1(idx)]';
    % %         ISImat_tst1=cell2mat(ISImat_tst1);
    % %         %tst2
    % %         ISImat_tst2=[sum([No_ISI1n(idx,:); No_ISI2n(idx,:); No_ISI3n(idx,:)]) sum([No_ISI4n(idx,:); No_ISI5n(idx,:); No_ISI6n(idx,:)]) No_ISI_tst2(idx)]';
    % %         ISImat_tst2=cell2mat(ISImat_tst2);
    % %         %tst3
    % %         ISImat_tst3=[sum([No_ISI1n(idx,:); No_ISI2n(idx,:); No_ISI3n(idx,:)]) sum([No_ISI4n(idx,:); No_ISI5n(idx,:); No_ISI6n(idx,:)]) No_ISI_tst3(idx)]';
    % %         ISImat_tst3=cell2mat(ISImat_tst3);
    % %         %tst4
    % %         ISImat_tst4=[sum([No_ISI1n(idx,:); No_ISI2n(idx,:); No_ISI3n(idx,:)]) sum([No_ISI4n(idx,:); No_ISI5n(idx,:); No_ISI6n(idx,:)]) No_ISI_tst4(idx)]';
    % %         ISImat_tst4=cell2mat(ISImat_tst4);
    % %         %tst5
    % %         ISImat_tst5=[sum([No_ISI1n(idx,:); No_ISI2n(idx,:); No_ISI3n(idx,:)]) sum([No_ISI4n(idx,:); No_ISI5n(idx,:); No_ISI6n(idx,:)]) No_ISI_tst5(idx)]';
    % %         ISImat_tst5=cell2mat(ISImat_tst5);
    % %         %tst6
    % %         ISImat_tst6=[sum([No_ISI1n(idx,:); No_ISI2n(idx,:); No_ISI3n(idx,:)]) sum([No_ISI4n(idx,:); No_ISI5n(idx,:); No_ISI6n(idx,:)]) No_ISI_tst6(idx)]';
    % %         ISImat_tst6=cell2mat(ISImat_tst6);
    % %tst1
    %         ISImat_tst1=[No_ISI1n(idx,:)  No_ISI4n(idx,:) No_ISI_tst1(idx)]';
    %         ISImat_tst1=cell2mat(ISImat_tst1);
    %         %tst2
    %         ISImat_tst2=[No_ISI1n(idx,:)  No_ISI4n(idx,:) No_ISI_tst2(idx)]';
    %         ISImat_tst2=cell2mat(ISImat_tst2);
    %         %tst3
    %         ISImat_tst3=[No_ISI1n(idx,:)  No_ISI4n(idx,:) No_ISI_tst3(idx)]';
    %         ISImat_tst3=cell2mat(ISImat_tst3);
    %         %tst4
    %         ISImat_tst4=[No_ISI1n(idx,:)  No_ISI4n(idx,:) No_ISI_tst4(idx)]';
    %         ISImat_tst4=cell2mat(ISImat_tst4);
    %         %tst5
    %         ISImat_tst5=[No_ISI1n(idx,:)  No_ISI4n(idx,:) No_ISI_tst5(idx)]';
    %         ISImat_tst5=cell2mat(ISImat_tst5);
    %         %tst6
    %         ISImat_tst6=[No_ISI1n(idx,:)  No_ISI4n(idx,:) No_ISI_tst6(idx)]';
    %         ISImat_tst6=cell2mat(ISImat_tst6);
    %         MLM=corr(ISImat_tst1')*length(find(ISImat_tst1>0));%Most Likely Motion
    %         MLMtst_1(idx,:)=MLM(3,:);
    %
    %         MLM=corr(ISImat_tst2')*length(find(ISImat_tst2>0));%Most Likely Motion
    %         MLMtst_2(idx,:)=MLM(3,:);
    %
    %         MLM=corr(ISImat_tst3')*length(find(ISImat_tst3>0));%Most Likely Motion
    %         MLMtst_3(idx,:)=MLM(3,:);
    %
    %         MLM=corr(ISImat_tst4')*length(find(ISImat_tst4>0));%Most Likely Motion
    %         MLMtst_4(idx,:)=MLM(3,:);
    %
    %         MLM=corr(ISImat_tst5')*length(find(ISImat_tst5>0));%Most Likely Motion
    %         MLMtst_5(idx,:)=MLM(3,:);
    %
    %         MLM=corr(ISImat_tst6')*length(find(ISImat_tst6>0));%Most Likely Motion
    %         MLMtst_6(idx,:)=MLM(3,:);
    %  end
    %     clear('MLM')
    %     RESULTS1(idx_ZTST,:)=((nansum(MLMtst_1)))/max((nansum(MLMtst_1)));
    %    RESULTS2(idx_ZTST,:)=((nansum(MLMtst_2)))/max((nansum(MLMtst_2)));
    %      RESULTS3(idx_ZTST,:)=((nansum(MLMtst_3)))/max((nansum(MLMtst_3)));
    %      RESULTS4(idx_ZTST,:)=((nansum(MLMtst_4)))/max((nansum(MLMtst_4)));
    %      RESULTS5(idx_ZTST,:)=((nansum(MLMtst_5)))/max((nansum(MLMtst_5)));
    %      RESULTS6(idx_ZTST,:)=((nansum(MLMtst_6)))/max((nansum(MLMtst_6)));
    for idx=1:length(idx_neuron);
        %         %tst1
        %         ISImat_tst1=[sum([No_ISI1n(idx,:); No_ISI2n(idx,:); No_ISI3n(idx,:)]) sum([No_ISI4n(idx,:); No_ISI5n(idx,:); No_ISI6n(idx,:)]) No_ISI_tst1(idx)]';
        %         ISImat_tst1=cell2mat(ISImat_tst1);
        %         %tst2
        %         ISImat_tst2=[sum([No_ISI1n(idx,:); No_ISI2n(idx,:); No_ISI3n(idx,:)]) sum([No_ISI4n(idx,:); No_ISI5n(idx,:); No_ISI6n(idx,:)]) No_ISI_tst2(idx)]';
        %         ISImat_tst2=cell2mat(ISImat_tst2);
        %         %tst3
        %         ISImat_tst3=[sum([No_ISI1n(idx,:); No_ISI2n(idx,:); No_ISI3n(idx,:)]) sum([No_ISI4n(idx,:); No_ISI5n(idx,:); No_ISI6n(idx,:)]) No_ISI_tst3(idx)]';
        %         ISImat_tst3=cell2mat(ISImat_tst3);
        %         %tst4
        %         ISImat_tst4=[sum([No_ISI1n(idx,:); No_ISI2n(idx,:); No_ISI3n(idx,:)]) sum([No_ISI4n(idx,:); No_ISI5n(idx,:); No_ISI6n(idx,:)]) No_ISI_tst4(idx)]';
        %         ISImat_tst4=cell2mat(ISImat_tst4);
        %         %tst5
        %         ISImat_tst5=[sum([No_ISI1n(idx,:); No_ISI2n(idx,:); No_ISI3n(idx,:)]) sum([No_ISI4n(idx,:); No_ISI5n(idx,:); No_ISI6n(idx,:)]) No_ISI_tst5(idx)]';
        %         ISImat_tst5=cell2mat(ISImat_tst5);
        %         %tst6
        %         ISImat_tst6=[sum([No_ISI1n(idx,:); No_ISI2n(idx,:); No_ISI3n(idx,:)]) sum([No_ISI4n(idx,:); No_ISI5n(idx,:); No_ISI6n(idx,:)]) No_ISI_tst6(idx)]';
        %         ISImat_tst6=cell2mat(ISImat_tst6);
        %tst1
        ISImat_tst1=[No_ISI_tst1(idx) ; zeros(1,135) ;No_ISI1n(idx,:);zeros(1,135);  No_ISI2n(idx,:) ;zeros(1,135);...
            No_ISI3n(idx,:);zeros(1,135);  No_ISI4n(idx,:) ;zeros(1,135);...
            No_ISI5n(idx,:);zeros(1,135);  No_ISI6n(idx,:) ;zeros(1,135)]';
        ISImat_tst1=cell2mat(ISImat_tst1);
        ISImat_tst1(isnan(ISImat_tst1))=0;
        %tst2
        ISImat_tst2=[ No_ISI_tst2(idx) ; zeros(1,135) ;No_ISI1n(idx,:);zeros(1,135);  No_ISI2n(idx,:) ;zeros(1,135);...
            No_ISI3n(idx,:);zeros(1,135);  No_ISI4n(idx,:) ;zeros(1,135);...
            No_ISI5n(idx,:);zeros(1,135);  No_ISI6n(idx,:) ;zeros(1,135)]';
        ISImat_tst2=cell2mat(ISImat_tst2);
        ISImat_tst2(isnan(ISImat_tst2))=0;
        %tst3
        ISImat_tst3=[No_ISI_tst3(idx); zeros(1,135) ;No_ISI1n(idx,:);zeros(1,135);  No_ISI2n(idx,:) ;zeros(1,135);...
            No_ISI3n(idx,:);zeros(1,135);  No_ISI4n(idx,:) ;zeros(1,135);...
            No_ISI5n(idx,:);zeros(1,135);  No_ISI6n(idx,:) ;zeros(1,135)]';
        ISImat_tst3=cell2mat(ISImat_tst3);
        ISImat_tst3(isnan(ISImat_tst3))=0;
        %tst4
        ISImat_tst4=[No_ISI_tst4(idx); zeros(1,135) ;No_ISI1n(idx,:);zeros(1,135);  No_ISI2n(idx,:) ;zeros(1,135);...
            No_ISI3n(idx,:);zeros(1,135);  No_ISI4n(idx,:) ;zeros(1,135);...
            No_ISI5n(idx,:);zeros(1,135);  No_ISI6n(idx,:) ;zeros(1,135)]';
        ISImat_tst4=cell2mat(ISImat_tst4);
        ISImat_tst4(isnan(ISImat_tst4))=0;
        %tst5
        ISImat_tst5=[No_ISI_tst5(idx); zeros(1,135) ;No_ISI1n(idx,:);zeros(1,135);  No_ISI2n(idx,:) ;zeros(1,135);...
            No_ISI3n(idx,:);zeros(1,135);  No_ISI4n(idx,:) ;zeros(1,135);...
            No_ISI5n(idx,:);zeros(1,135);  No_ISI6n(idx,:) ;zeros(1,135)]';
        ISImat_tst5=cell2mat(ISImat_tst5);
        
        ISImat_tst5(isnan(ISImat_tst5))=0;
        %tst6
        ISImat_tst6=[No_ISI_tst6(idx); zeros(1,135) ;No_ISI1n(idx,:);zeros(1,135);  No_ISI2n(idx,:) ;zeros(1,135);...
            No_ISI3n(idx,:);zeros(1,135);  No_ISI4n(idx,:) ;zeros(1,135);...
            No_ISI5n(idx,:);zeros(1,135);  No_ISI6n(idx,:) ;zeros(1,135)]';
        ISImat_tst6=cell2mat(ISImat_tst6);
        ISImat_tst6(isnan(ISImat_tst6))=0;
        %auto correlation using cross correlation with the test data with
        %itself and training data. a=correlation from 0 to 1 b = lag
        %tst1 autocorr
        [a b]=xcorr(interp(cell2mat(No_ISI_tst1(:,idx)),3),interp(ISImat_tst1,3));
        if max(a)>0; %ignore non firing neurons
            [pks idx_pks]=findpeaks(a/max(a),'MinPeakDistance',2*135);
            idx_pks_itself=find(idx_pks<(135*43)&idx_pks>(135*41));
            idx_pks_raddev=find(idx_pks<(135*37)&idx_pks>(135*32));
            idx_pks_wstflx=find(idx_pks<(135*19)&idx_pks>(135*14));
            idx_pks_raddev1kg=find(idx_pks<(135*31)&idx_pks>(135*26));
            idx_pks_wstflx1kg=find(idx_pks<(135*13)&idx_pks>(135*8));
            idx_pks_raddev2kg=find(idx_pks<(135*25)&idx_pks>(135*20));
            idx_pks_wstflx2kg=find(idx_pks<(135*7)&idx_pks>(135*2));
            cr_raddev=max(pks(idx_pks_raddev));
            cr_wstflx=max(pks(idx_pks_wstflx));
            cr_raddev1kg=max(pks(idx_pks_raddev1kg));
            cr_wstflx1kg=max(pks(idx_pks_wstflx1kg));
            cr_raddev2kg=max(pks(idx_pks_raddev2kg));
            cr_wstflx2kg=max(pks(idx_pks_wstflx2kg));
            cr_itself=max(pks(idx_pks_itself));
            
            MLMtst_1(idx,:)=[cr_raddev cr_raddev1kg cr_raddev2kg cr_wstflx cr_wstflx1kg cr_wstflx2kg ];
        end
        clear('a','b','pks','idx_pks','idx_pks_itself','idx_pks_raddev','idx_pks_wstflx');
        %tst2 autocorr
        [a b]=xcorr(interp(cell2mat(No_ISI_tst2(:,idx)),3),interp(ISImat_tst2,3));
        if max(a)>0;
            [pks idx_pks]=findpeaks(a/max(a),'MinPeakDistance',2*135);
            idx_pks_itself=find(idx_pks<(135*43)&idx_pks>(135*41));
            idx_pks_raddev=find(idx_pks<(135*37)&idx_pks>(135*32));
            idx_pks_wstflx=find(idx_pks<(135*19)&idx_pks>(135*14));
            idx_pks_raddev1kg=find(idx_pks<(135*31)&idx_pks>(135*26));
            idx_pks_wstflx1kg=find(idx_pks<(135*13)&idx_pks>(135*8));
            idx_pks_raddev2kg=find(idx_pks<(135*25)&idx_pks>(135*20));
            idx_pks_wstflx2kg=find(idx_pks<(135*7)&idx_pks>(135*2));
            cr_raddev=max(pks(idx_pks_raddev));
            cr_wstflx=max(pks(idx_pks_wstflx));
            cr_raddev1kg=max(pks(idx_pks_raddev1kg));
            cr_wstflx1kg=max(pks(idx_pks_wstflx1kg));
            cr_raddev2kg=max(pks(idx_pks_raddev2kg));
            cr_wstflx2kg=max(pks(idx_pks_wstflx2kg));
            cr_itself=max(pks(idx_pks_itself));
            MLMtst_2(idx,:)=[cr_raddev cr_raddev1kg cr_raddev2kg cr_wstflx cr_wstflx1kg cr_wstflx2kg ];
        end
        clear('a','b','pks','idx_pks','idx_pks_itself','idx_pks_raddev','idx_pks_wstflx');
        
        %tst3 autocorr
        [a b]=xcorr(interp(cell2mat(No_ISI_tst3(:,idx)),3),interp(ISImat_tst3,3));
        if max(a)>0;
            [pks idx_pks]=findpeaks(a/max(a),'MinPeakDistance',2*135);
            idx_pks_itself=find(idx_pks<(135*43)&idx_pks>(135*41));
            idx_pks_raddev=find(idx_pks<(135*37)&idx_pks>(135*32));
            idx_pks_wstflx=find(idx_pks<(135*19)&idx_pks>(135*14));
            idx_pks_raddev1kg=find(idx_pks<(135*31)&idx_pks>(135*26));
            idx_pks_wstflx1kg=find(idx_pks<(135*13)&idx_pks>(135*8));
            idx_pks_raddev2kg=find(idx_pks<(135*25)&idx_pks>(135*20));
            idx_pks_wstflx2kg=find(idx_pks<(135*7)&idx_pks>(135*2));
            cr_raddev=max(pks(idx_pks_raddev));
            cr_wstflx=max(pks(idx_pks_wstflx));
            cr_raddev1kg=max(pks(idx_pks_raddev1kg));
            cr_wstflx1kg=max(pks(idx_pks_wstflx1kg));
            cr_raddev2kg=max(pks(idx_pks_raddev2kg));
            cr_wstflx2kg=max(pks(idx_pks_wstflx2kg));
            cr_itself=max(pks(idx_pks_itself));
            MLMtst_3(idx,:)=[cr_raddev cr_raddev1kg cr_raddev2kg cr_wstflx cr_wstflx1kg cr_wstflx2kg ];
        end
        clear('a','b','pks','idx_pks','idx_pks_itself','idx_pks_raddev','idx_pks_wstflx');
        %tst4 autocorr
        [a b]=xcorr(interp(cell2mat(No_ISI_tst4(:,idx)),3),interp(ISImat_tst4,3));
        if max(a)>0;
            [pks idx_pks]=findpeaks(a/max(a),'MinPeakDistance',2*135);
            idx_pks_itself=find(idx_pks<(135*43)&idx_pks>(135*41));
            idx_pks_raddev=find(idx_pks<(135*37)&idx_pks>(135*32));
            idx_pks_wstflx=find(idx_pks<(135*19)&idx_pks>(135*14));
            idx_pks_raddev1kg=find(idx_pks<(135*31)&idx_pks>(135*26));
            idx_pks_wstflx1kg=find(idx_pks<(135*13)&idx_pks>(135*8));
            idx_pks_raddev2kg=find(idx_pks<(135*25)&idx_pks>(135*20));
            idx_pks_wstflx2kg=find(idx_pks<(135*7)&idx_pks>(135*2));
            cr_raddev=max(pks(idx_pks_raddev));
            cr_wstflx=max(pks(idx_pks_wstflx));
            cr_raddev1kg=max(pks(idx_pks_raddev1kg));
            cr_wstflx1kg=max(pks(idx_pks_wstflx1kg));
            cr_raddev2kg=max(pks(idx_pks_raddev2kg));
            cr_wstflx2kg=max(pks(idx_pks_wstflx2kg));
            cr_itself=max(pks(idx_pks_itself));
            MLMtst_4(idx,:)=[cr_raddev cr_raddev1kg cr_raddev2kg cr_wstflx cr_wstflx1kg cr_wstflx2kg ];
        end
        clear('a','b','pks','idx_pks','idx_pks_itself','idx_pks_raddev','idx_pks_wstflx');
        %tst5 autocorr
        [a b]=xcorr(interp(cell2mat(No_ISI_tst5(:,idx)),3),interp(ISImat_tst5,3));
        if max(a)>0;
            [pks idx_pks]=findpeaks(a/max(a),'MinPeakDistance',2*135);
            idx_pks_itself=find(idx_pks<(135*43)&idx_pks>(135*41));
            idx_pks_raddev=find(idx_pks<(135*37)&idx_pks>(135*32));
            idx_pks_wstflx=find(idx_pks<(135*19)&idx_pks>(135*14));
            idx_pks_raddev1kg=find(idx_pks<(135*31)&idx_pks>(135*26));
            idx_pks_wstflx1kg=find(idx_pks<(135*13)&idx_pks>(135*8));
            idx_pks_raddev2kg=find(idx_pks<(135*25)&idx_pks>(135*20));
            idx_pks_wstflx2kg=find(idx_pks<(135*7)&idx_pks>(135*2));
            cr_raddev=max(pks(idx_pks_raddev));
            cr_wstflx=max(pks(idx_pks_wstflx));
            cr_raddev1kg=max(pks(idx_pks_raddev1kg));
            cr_wstflx1kg=max(pks(idx_pks_wstflx1kg));
            cr_raddev2kg=max(pks(idx_pks_raddev2kg));
            cr_wstflx2kg=max(pks(idx_pks_wstflx2kg));
            cr_itself=max(pks(idx_pks_itself));
            MLMtst_5(idx,:)=[cr_raddev cr_raddev1kg cr_raddev2kg cr_wstflx cr_wstflx1kg cr_wstflx2kg ];
        end
        clear('a','b','pks','idx_pks','idx_pks_itself','idx_pks_raddev','idx_pks_wstflx');
        %tst6 autocorr
        [a b]=xcorr(interp(cell2mat(No_ISI_tst6(:,idx)),3),interp(ISImat_tst6,3));
        if max(a)>0;
            [pks idx_pks]=findpeaks(a/max(a),'MinPeakDistance',2*135);
            idx_pks_itself=find(idx_pks<(135*43)&idx_pks>(135*41));
            idx_pks_raddev=find(idx_pks<(135*37)&idx_pks>(135*32));
            idx_pks_wstflx=find(idx_pks<(135*19)&idx_pks>(135*14));
            idx_pks_raddev1kg=find(idx_pks<(135*31)&idx_pks>(135*26));
            idx_pks_wstflx1kg=find(idx_pks<(135*13)&idx_pks>(135*8));
            idx_pks_raddev2kg=find(idx_pks<(135*25)&idx_pks>(135*20));
            idx_pks_wstflx2kg=find(idx_pks<(135*7)&idx_pks>(135*2));
            cr_raddev=max(pks(idx_pks_raddev));
            cr_wstflx=max(pks(idx_pks_wstflx));
            cr_raddev1kg=max(pks(idx_pks_raddev1kg));
            cr_wstflx1kg=max(pks(idx_pks_wstflx1kg));
            cr_raddev2kg=max(pks(idx_pks_raddev2kg));
            cr_wstflx2kg=max(pks(idx_pks_wstflx2kg));
            cr_itself=max(pks(idx_pks_itself));
            MLMtst_6(idx,:)=[cr_raddev cr_raddev1kg cr_raddev2kg cr_wstflx cr_wstflx1kg cr_wstflx2kg ];
        end
        clear('a','b','pks','idx_pks','idx_pks_itself','idx_pks_raddev','idx_pks_wstflx');
    end
    clear('MLM')
    %clear non firing neurons
    MLMtst_1([find(MLMtst_1(:,1)==0 & MLMtst_1(:,2)==0 & MLMtst_1(:,3)==0 & MLMtst_1(:,4)==0 & MLMtst_1(:,6)==0 & MLMtst_1(:,6)==0)],:)=[];
    MLMtst_2([find(MLMtst_2(:,1)==0 & MLMtst_2(:,2)==0 & MLMtst_2(:,3)==0 & MLMtst_2(:,4)==0 & MLMtst_2(:,6)==0 & MLMtst_2(:,6)==0)],:)=[];
    MLMtst_3([find(MLMtst_3(:,1)==0 & MLMtst_3(:,2)==0 & MLMtst_3(:,3)==0 & MLMtst_3(:,4)==0 & MLMtst_3(:,6)==0 & MLMtst_3(:,6)==0)],:)=[];
    MLMtst_4([find(MLMtst_4(:,1)==0 & MLMtst_4(:,2)==0 & MLMtst_4(:,3)==0 & MLMtst_4(:,4)==0 & MLMtst_4(:,6)==0 & MLMtst_4(:,6)==0)],:)=[];
    MLMtst_5([find(MLMtst_5(:,1)==0 & MLMtst_5(:,2)==0 & MLMtst_5(:,3)==0 & MLMtst_5(:,4)==0 & MLMtst_5(:,6)==0 & MLMtst_5(:,6)==0)],:)=[];
    MLMtst_6([find(MLMtst_6(:,1)==0 & MLMtst_6(:,2)==0 & MLMtst_6(:,3)==0 & MLMtst_6(:,4)==0 & MLMtst_6(:,6)==0 & MLMtst_6(:,6)==0)],:)=[];
    
    
    RESULTS1(idx_ZTST,:)=median(MLMtst_1);
    RESULTS2(idx_ZTST,:)=median(MLMtst_2);
    RESULTS3(idx_ZTST,:)=median(MLMtst_3);
    RESULTS4(idx_ZTST,:)=median(MLMtst_4);
    RESULTS5(idx_ZTST,:)=median(MLMtst_5);
    RESULTS6(idx_ZTST,:)=median(MLMtst_6);
    
    clear('MLMtst_6', 'MLMtst_1', 'MLMtst_2', 'MLMtst_3', 'MLMtst_4', 'MLMtst_5');
end
%% Results and plots

% ENG 
figure(40)
plot((1:length(all_data))/fs_all,all_data(102,:))
ylabel('Nerve ENG (Volts)');
xlabel('Time(seconds)');
title('ENG with bandpass and upsampled');
grid minor;

% Muscle activity per motion over time under 3 weights

% Radial Deviation (ECR)
figure(31)
plot((1:wndw)/fs_all,smooth(RadialDeviation(idx_EXTCRAD_EMG,:),400),'r'); hold on;
plot((wndw+1:2*wndw)/fs_all,smooth(RadialDeviation_1kg(idx_EXTCRAD_EMG,:),400),'g'); hold on;
plot((2*wndw+1:3*wndw)/fs_all,smooth(RadialDeviation_2kg(idx_EXTCRAD_EMG,:),400),'b');
ylabel('RMS Muscle Activity');
xlabel('Time(seconds)');
title('RadialDeviation then WristFlexion');
hold on;

% Wrist Flexion (FCR)
plot((3*wndw+1:4*wndw)/fs_all,smooth(WristFlexion(idx_FLXCRAD_EMG,:),400),'r');hold on;
plot((4*wndw+1:5*wndw)/fs_all,smooth(WristFlexion_1kg(idx_FLXCRAD_EMG,:),400),'g');hold on;
plot((5*wndw+1:6*wndw)/fs_all,smooth(WristFlexion_2kg(idx_FLXCRAD_EMG,:),400),'b');hold on;
hold on;
legend('ECR','ECR1KG','ECR2KG','FCR','FCR1KG','FCR2KG');
set(legend,'Location','northeast');
grid minor;

%IMU data
figure(32)
for i=1:3
    hold on;
    plot((1:length(all_data_DORSUM_Acc))/fs_all,smooth(all_data_DORSUM_Acc(i,:)*((length(Templates))*-1)+3,1200)); hold on
end
ylabel('Arbitrary Value of Deviation');
xlabel('Time(seconds)');
title('RadialDeviation then WristFlexion');
hold on;
legend('X','Y','Z');

% All templates over time
clour=char('k','m', 'c' ,'r', 'g', 'b', 'w', 'y' ,'r--', 'b--', 'g--');
figure(33)
for idx=1:length(Templates(:,1))
    % Collect indices per spike label
    idx_lbl{idx,:}=find(MUAP_label==idx);
    for i=1:length(Templates{idx,:}(:,1));
        % Create cell for true indices per spike label
        idx_lbl{idx,:}(:,i)=MUAP_idx(idx_lbl{idx,:}(:,i));
        % Plot all templates
        wndow=((idx_lbl{idx,:}(:,i))+(1:length(Templates{idx,:}(i,:)))-1);
        plot(wndow/fs_all,idx+Templates{idx,:}(i,:),clour(idx)); hold on; %Each row corresponds to a spike label
    end
end
clear('i', 'idx_eur', 'wndow');
hold on;
plot((1:length(all_data))/fs_all, scrs_all(2,:)*length(Templates),'k','linewidth',1.25); % basic imu+emg
ylabel('Spike Coloured per Neuron');
xlabel('Time(seconds)');
ylim([-3 16]);
title('RadialDeviation then WristFlexion');
hold on;
grid minor;

%ISI dist examples
figure(34)
plot((1:3*length(No_ISI1n))/3000,abs(interp(((nansum(No_ISI1n))/max(nansum(No_ISI1n))),3)),'b','linewidth',1.25); hold on
plot((1:3*length(No_ISI1n))/3000,abs(interp(((nansum(cell2mat(No_ISI_tst1')))/max(nansum(cell2mat(No_ISI_tst1')))),3)),'r','linewidth',1)
ylim([0 1]);
xlim([0.005 0.135]);
ylabel('Firing Probability');
xlabel('Interspike Interval (Seconds)');
title('Radial Deviation 0kg ISI Distribution');
legend('Train','Test');
set(legend,'Location','northeast');
grid minor
figure(35)
plot((1:3*length(No_ISI1n))/3000,abs(interp(((nansum(No_ISI1n))/max(nansum(No_ISI1n))),3)),'b','linewidth',1.25); hold on
plot((1:3*length(No_ISI1n))/3000,abs(interp(((nansum(cell2mat(No_ISI_tst2')))/max(nansum(cell2mat(No_ISI_tst2')))),3)),'r','linewidth',1)
ylim([0 1]);
xlim([0.005 0.135]);
ylabel('Firing Probability');
xlabel('Interspike Interval (Seconds)');
title('Radial Deviation 1kg ISI Distribution');
legend('Train','Test');
set(legend,'Location','northeast');
grid minor
figure(36)
plot((1:3*length(No_ISI1n))/3000,abs(interp(((nansum(No_ISI1n))/max(nansum(No_ISI1n))),3)),'b','linewidth',1.25); hold on
plot((1:3*length(No_ISI1n))/3000,abs(interp(((nansum(cell2mat(No_ISI_tst3')))/max(nansum(cell2mat(No_ISI_tst3')))),3)),'r','linewidth',1)
ylim([0 1]);
xlim([0.005 0.135]);
ylabel('Firing Probability');
xlabel('Interspike Interval (Seconds)');
title('Radial Deviation 2kg ISI Distribution');
legend('Train','Test');
set(legend,'Location','northeast');
grid minor
figure(37)
plot((1:3*length(No_ISI4n))/3000,abs(interp(((nansum(No_ISI4n))/max(nansum(No_ISI4n))),3)),'b','linewidth',1.25); hold on
plot((1:3*length(No_ISI4n))/3000,abs(interp(((nansum(cell2mat(No_ISI_tst4')))/max(nansum(cell2mat(No_ISI_tst4')))),3)),'r','linewidth',1)
ylim([0 1]);
xlim([0.005 0.135]);
ylabel('Firing Probability');
xlabel('Interspike Interval (Seconds)');
title('Wrist Flexion 0kg ISI Distribution');
legend('Train','Test');
set(legend,'Location','northeast');
grid minor
figure(38)
plot((1:3*length(No_ISI4n))/3000,abs(interp(((nansum(No_ISI4n))/max(nansum(No_ISI4n))),3)),'b','linewidth',1.25); hold on
plot((1:3*length(No_ISI4n))/3000,abs(interp(((nansum(cell2mat(No_ISI_tst5')))/max(nansum(cell2mat(No_ISI_tst5')))),3)),'r','linewidth',1)
ylim([0 1]);
xlim([0.005 0.135]);
ylabel('Firing Probability');
xlabel('Interspike Interval (Seconds)');
title('Wrist Flexion 1kg ISI Distribution');
legend('Train','Test');
set(legend,'Location','northeast');
grid minor
figure(39)
plot((1:3*length(No_ISI4n))/3000,abs(interp(((nansum(No_ISI4n))/max(nansum(No_ISI4n))),3)),'b','linewidth',1.25); hold on
plot((1:3*length(No_ISI4n))/3000,abs(interp(((nansum(cell2mat(No_ISI_tst6')))/max(nansum(cell2mat(No_ISI_tst6')))),3)),'r','linewidth',1)
ylim([0 1]);
xlim([0.005 0.135]);
ylabel('Firing Probability');
xlabel('Interspike Interval (Seconds)');
title('Wrist Flexion 2kg ISI Distribution');
legend('Train','Test');
set(legend,'Location','northeast');
grid minor

%Most likely motion (Correlation with median, maximum, minimum and
%quantiles)
figure(21);
boxplot(RESULTS1); hold on;
ylabel('Correlation'); xlim([0.5 6.5]); ylim([0 1]);
xlabel('Motion Index (1=Rad0kg 2=Rad1kg 3=Rad2kg 4=Wrist0kg 5=Wrist1kg 6=Wrist2kg)');
title('Motion 1: Radial Deviation 0Kg');
grid('minor');
figure(22);
boxplot(RESULTS2); hold on;
ylabel('Correlation'); xlim([0.5 6.5]); ylim([0 1]);
xlabel('Motion Index (1=Rad0kg 2=Rad1kg 3=Rad2kg 4=Wrist0kg 5=Wrist1kg 6=Wrist2kg)');
title('Motion 2: Radial Deviation 1Kg');
grid('minor');
figure(23);
boxplot(RESULTS3); hold on;
ylabel('Correlation'); xlim([0.5 6.5]); ylim([0 1]);
xlabel('Motion Index (1=Rad0kg 2=Rad1kg 3=Rad2kg 4=Wrist0kg 5=Wrist1kg 6=Wrist2kg)');
title('Motion 3: Radial Deviation 2Kg');
grid('minor');
figure(24);
boxplot(RESULTS4); hold on;
ylabel('Correlation'); xlim([0.5 6.5]); ylim([0 1]);
xlabel('Motion Index (1=Rad0kg 2=Rad1kg 3=Rad2kg 4=Wrist0kg 5=Wrist1kg 6=Wrist2kg)');
title('Motion 4: Wrist Flexion 0Kg');
grid('minor');
figure(25);
boxplot(RESULTS5); hold on;
ylabel('Correlation'); xlim([0.5 6.5]); ylim([0 1]);
xlabel('Motion Index (1=Rad0kg 2=Rad1kg 3=Rad2kg 4=Wrist0kg 5=Wrist1kg 6=Wrist2kg)');
title('Motion 5: Wrist Flexion 1Kg');
grid('minor');
figure(26);
boxplot(RESULTS6); hold on;
ylabel('Correlation'); xlim([0.5 6.5]); ylim([0 1]);
xlabel('Motion Index (1=Rad0kg 2=Rad1kg 3=Rad2kg 4=Wrist0kg 5=Wrist1kg 6=Wrist2kg)');
title('Motion 6: Wrist Flexion 2Kg');
grid('minor');


%Precision and recall
% tp=sum(sum([RESULTS1(:,1) RESULTS2(:,2) RESULTS3(:,3) RESULTS4(:,4) RESULTS5(:,5) RESULTS6(:,6)]))/30
% RESULTS1(:,1)=[];
% RESULTS2(:,2)=[] ;
% RESULTS3(:,3)=[] ;
% RESULTS4(:,4)=[] ;
% RESULTS5(:,5)=[] ;
% RESULTS6(:,6)=[];
% RESULTS1(:,6)=[];
% RESULTS2(:,6)=[] ;
% RESULTS3(:,6)=[] ;
% RESULTS4(:,6)=[] ;
% RESULTS5(:,6)=[] ;
% RESULTS6(:,6)=[];
% fp=sum(sum(abs([RESULTS1 RESULTS2 RESULTS3 RESULTS4 RESULTS5 RESULTS6])))/150
% precsn=tp/(fp+tp)
% tp=sum(sum([RESULTS1(:,1) RESULTS2(:,1) RESULTS3(:,1) RESULTS4(:,2) RESULTS5(:,2) RESULTS6(:,2)]))/30
% RESULTS1(:,1)=[];
% RESULTS2(:,1)=[] ;
% RESULTS3(:,1)=[] ;
% RESULTS4(:,2)=[] ;
% RESULTS5(:,2)=[] ;
% RESULTS6(:,2)=[];
% % RESULTS1(:,2)=[];
% % RESULTS2(:,2)=[] ;
% % RESULTS3(:,2)=[] ;
% % RESULTS4(:,2)=[] ;
% % RESULTS5(:,2)=[] ;
% % RESULTS6(:,2)=[];
% fp=sum(sum([RESULTS1 RESULTS2 RESULTS3 RESULTS4 RESULTS5 RESULTS6]))/30
% precsn=tp/(fp+tp)
tp=sum(sum([RESULTS1(:,1) RESULTS2(:,2) RESULTS3(:,3) RESULTS4(:,4) RESULTS5(:,5) RESULTS6(:,6)]))/30
RESULTS1(:,1)=[];
RESULTS2(:,2)=[] ;
RESULTS3(:,3)=[] ;
RESULTS4(:,4)=[] ;
RESULTS5(:,5)=[] ;
% RESULTS6(:,6)=[];
% RESULTS1(:,6)=[];
% RESULTS2(:,6)=[] ;
% RESULTS3(:,6)=[] ;
% RESULTS4(:,6)=[] ;
% RESULTS5(:,6)=[] ;
% RESULTS6(:,6)=[];
fp=sum(sum(abs([RESULTS1 RESULTS2 RESULTS3 RESULTS4 RESULTS5 RESULTS6])))/150
precsn=tp/(fp+tp)

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

%% Classify directedness/smoothness

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
