%% Cumulative muscle activity per motion over time under 3 weights 

% Radial Deviation
figure()
plot(cumsum(RadialDeviation(idx_FLXCRAD_EMG,:)),'r--');hold on;
plot(cumsum(RadialDeviation(idx_EXTCRAD_EMG,:)),'b--');
plot(cumsum(RadialDeviation(idx_EXTCULN_EMG,:)),'g--');
plot(cumsum(RadialDeviation_1kg(idx_FLXCRAD_EMG,:)),'r-.');hold on;
plot(cumsum(RadialDeviation_1kg(idx_EXTCRAD_EMG,:)),'b-.');
plot(cumsum(RadialDeviation_1kg(idx_EXTCULN_EMG,:)),'g-.');
plot(cumsum(RadialDeviation_2kg(idx_FLXCRAD_EMG,:)),'r-');hold on;
plot(cumsum(RadialDeviation_2kg(idx_EXTCRAD_EMG,:)),'b-');
plot(cumsum(RadialDeviation_2kg(idx_EXTCULN_EMG,:)),'g-');
ylabel('Cumulative Muscle Activity');
xlabel('Time(samples)');
title('RadialDeviation');
hold on;
legend('FCR','ECR','ECU','FCR1KG','ECR1KG','ECU1KG','FCR2KG','ECR2KG','ECU2KG');
set(legend,'FontSize',13,'Location','northwest');

% Wrist Flexion
figure()
plot(cumsum(WristFlexion(idx_FLXCRAD_EMG,:)),'r--');hold on;
plot(cumsum(WristFlexion(idx_EXTCRAD_EMG,:)),'b--');
plot(cumsum(WristFlexion(idx_EXTCULN_EMG,:)),'g--');
plot(cumsum(WristFlexion_1kg(idx_FLXCRAD_EMG,:)),'r-.');hold on;
plot(cumsum(WristFlexion_1kg(idx_EXTCRAD_EMG,:)),'b-.');
plot(cumsum(WristFlexion_1kg(idx_EXTCULN_EMG,:)),'g-.');
plot(cumsum(WristFlexion_2kg(idx_FLXCRAD_EMG,:)),'r-');hold on;
plot(cumsum(WristFlexion_2kg(idx_EXTCRAD_EMG,:)),'b-');
plot(cumsum(WristFlexion_2kg(idx_EXTCULN_EMG,:)),'g-');
ylabel('Cumulative Muscle Activity');
xlabel('Time(samples)');
title('WristFlexion');
hold on;
legend('FCR','ECR','ECU','FCR1KG','ECR1KG','ECU1KG','FCR2KG','ECR2KG','ECU2KG');
set(legend,'FontSize',13,'Location','northwest');

%% Spike Identiy vs Time and RMS muscle activity

scatter(MUAP_idx/fs_all,MUAP_label); hold on;
plot((1:length(all_RadialDeviation))/fs_all,((all_RadialDeviation(idx_EXTCRAD_EMG,...
    :))/max(all_RadialDeviation(idx_EXTCRAD_EMG,...
:)))*max(MUAP_label));
ylabel('SpikeID_Label');
xlabel('Time(seconds)');
title('WristFlexion'); %title('RadialDeviation');
grid minor;