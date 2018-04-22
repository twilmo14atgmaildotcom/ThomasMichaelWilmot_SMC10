%% Guitar Playing and Motor Interaction
% Thomas Michael Wilmot
% Student of Sound and Music Computing MSc at Aalborg Univeristy
% Last updated: 22/04/18
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
cd('C:\Users\Thomas\Documents\GitHub\SMC10_ThomasWilmot\SMC10ThomasWilmot\Delsysdata210318\TomSimpleMotionDelsys_210318');

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

% Remove old indices
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

%% Create labelled time windows from EMG and IMU data

%Detection method used for detecting spike onsets, offsets and duration did
%not translate for the following. Meaning labelled time windows could not
%be created that accurately represented the following without manually
%adjusting thresholds per sensor data and removing movement artifacts in
%data.

% Start of contraction (EMG)

% End of contraction (EMG)

% Rest position (minimum muscle activity)

% Held pose position (prolonged muscle activity)

% Onsets of motion (IMU)

% Offsets of motion (IMU)

% Percentage muscle strength used (3EMG)

% Weight used during motion (none, 1kg or 2kg)


%% Horizontally concatonate recordings with the same motions

all_RadialDeviation=horzcat(RadialDeviation, RadialDeviation_1kg, RadialDeviation_2kg);
all_WristFlexion=horzcat(WristFlexion, WristFlexion_1kg, WristFlexion_2kg);

%% Save individual nerve sensor data for same motions across all recordings

name_var='all_RadialDeviation';
curr_var=eval(name_var);
savefile=strcat(name_var,'_','NRVPOSTOUT_EMG.mat');
x=curr_var(idx_NRVPOSTOUT_EMG,:);
save(savefile,'x');
savefile=strcat(name_var,'_','NRVPOSTINT_EMG.mat');
x=curr_var(idx_NRVPOSTINT_EMG,:);
save(savefile,'x');
savefile=strcat(name_var,'_','NRVANTINT_EMG.mat');
x=curr_var(idx_NRVANTINT_EMG,:);
save(savefile,'x');
savefile=strcat(name_var,'_','NRVANTOUT_EMG.mat');
x=curr_var(idx_NRVANTOUT_EMG,:);
save(savefile,'x');


name_var='all_WristFlexion';
curr_var=eval(name_var);
savefile=strcat(name_var,'_','NRVPOSTOUT_EMG.mat');
x=curr_var(idx_NRVPOSTOUT_EMG,:);
save(savefile,'x');
savefile=strcat(name_var,'_','NRVPOSTINT_EMG.mat');
x=curr_var(idx_NRVPOSTINT_EMG,:);
save(savefile,'x');
savefile=strcat(name_var,'_','NRVANTINT_EMG.mat');
x=curr_var(idx_NRVANTINT_EMG,:);
save(savefile,'x');
savefile=strcat(name_var,'_','NRVANTOUT_EMG.mat');
x=curr_var(idx_NRVANTOUT_EMG,:);
save(savefile,'x');

clear('savefile', 'x', 'name_var', 'curr_var');

%% Save individual nerve sensor data

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

% run LEMG_Analyzer.m
% open desired folder
% enter sample frequency
% click on start clustering

%% Reference EMG data (if required)

%% Score spikes within time windows according to muscle and motion data

%% Compare spike counts from similar motions over increased weight applied during motion

%Repeat below for each motion
% C={[MUAP_label]};
% [uvals,~,bin]=unique(cell2mat(C));
% origarray=repelem(1:numel(C),cellfun(@numel,C));
% dist=table(uvals',accumarray(bin,1),accumarray(bin,origarray,[],@(x) {x}),'VariableNames',{'spikeID_label','repetition','arrayindices'});
% dist=sortrows(dist,'repetition');
