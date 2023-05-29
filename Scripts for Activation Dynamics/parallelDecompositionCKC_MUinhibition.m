    clear all
    close all
    clc 
    muscleChannels = 1:64; % 1:64 TA, soleus 65:128, GM - 129, GL - 130, PL - 131
    fs = 2048;
    exFactor = 16;
    LP = 20;
    HP = 500;
    NITER = 100; % 150
    signalQuality = 5;
    DiffModeChkbox = 0;
    iterPlots = false;
    costFunc = 5;
    plotFlag = 1;
    nRep = 8;
    color = '#A2142F';
    addpath(genpath('DEMUSE'))

    if muscleChannels(1) == 1
        mus = 'TA';
    elseif muscleChannels(1) == 65
        mus = 'SOL';
    elseif muscleChannels(1) == 129
        mus = 'GM';
    elseif muscleChannels(1) == 130
        mus = 'GL';
    else
        mus = 'PL';
    end

    %% Data extraction: MU inhibition (pilot)
   
    trial = 'Trial5C1R';
    trial1 = 'Trial3C1R1';
    path = 'C:\Users\User\Desktop\Thesis\Data\subject01\';
    MVC_Refa = 'A08_EMGMVC_Refa\';
% 
   EMG_MVC_1 = load(strcat(path,strcat(MVC_Refa,trial1)));             %load MVC trial
   EMG_MVC = EMG_MVC_1.OutputData.EMGData(muscleChannels,:);
  

   % Subject 6 has 3 MVC trials, so take average of the trials

%     EMG_MVC_1 = load('C:\Users\User\Desktop\Thesis\Data\subject06\A08_EMGMVC_Refa\Trial3C1R1.mat'); % load MVC trial1 
%     EMG_MVC_2 = load('C:\Users\User\Desktop\Thesis\Data\subject06\A08_EMGMVC_Refa\Trial3C1R2.mat'); % load MVC trial2
%     EMG_MVC_3 = load('C:\Users\User\Desktop\Thesis\Data\subject06\A08_EMGMVC_Refa\Trial3C1R3.mat'); % load MVC trial3
%     EMG_MVC_4 = EMG_MVC_1.OutputData.EMGData(muscleChannels,:);%1:length(EMG_MVC_6));
%     EMG_MVC_5 = zeros(length(muscleChannels),length(EMG_MVC_4));
%     EMG_MVC_6 = zeros(length(muscleChannels),length(EMG_MVC_4));
%     EMG_MVC_5(muscleChannels,1:length(EMG_MVC_2.OutputData.EMGData)) = EMG_MVC_2.OutputData.EMGData(muscleChannels,:);%1:length(EMG_MVC_6));
%     EMG_MVC_6(muscleChannels,1:length(EMG_MVC_3.OutputData.EMGData)) = EMG_MVC_3.OutputData.EMGData(muscleChannels,:);    
% 
%     for i = 1:64       %1 - bipolar, 64 - sol, TA
%         for j = 1:length(EMG_MVC_4)
%             EMG_MVC(i,j) = (EMG_MVC_4(i,j)+EMG_MVC_5(i,j)+EMG_MVC_6(i,j))/3; % MVC mean
%         end
%     end


    %% Data extraction

    EMG = struct('raw',[]);
   
for i = 1: nRep
    load(strcat(path,strcat(trial,num2str(i))));             %load actual trial
    % Synchronization
    % Parameters
    Achilles.rawAchData = OutputData.AchData;
    EMG(i).raw = OutputData.EMGData';                              % Raw EMG data
    Achilles.rawVar.pos = OutputData.AchData(:,2);              % Achilles raw position
    Achilles.rawVar.vel = OutputData.AchData(:,3);              % Achilles raw velocity
    Achilles.rawVar.tor = OutputData.AchData(:,4);              % Achilles raw torque
    Param.fs = 2048;                                            % Sampling frequency
    Param.nSamples = size(Achilles.rawAchData,1);               % Number of samples
    Achilles.rawVar.time = (0:1/Param.fs:(Param.nSamples-1)/Param.fs)';       % New time vector to avoid clock errors

    % Synchronization of Achilles and EMG data
    Param.trigAchOn = find(Achilles.rawAchData(:,8) == 1);      % Valid data from Achilles (only when fader factor value is 1)
    %Param.trigAchOn= ones(size(Achilles.rawAchData(:,8)))
    maxIndex = length(Achilles.rawAchData(:,1)) - Param.trigAchOn(end);
%     EMGsync  = syncEMGandAchilles(EMG(i), Achilles.rawVar.tor(Param.trigAchOn), maxIndex );
    EMGsync  = syncEMGandAchilles(EMG(i), Achilles.rawVar.tor(Param.trigAchOn));
    EMG(i).sync = EMGsync.sync;
    EMG(i).pos= Achilles.rawVar.pos(Param.trigAchOn);
    EMG(i).torque= Achilles.rawVar.tor(Param.trigAchOn);
    EMG(i).velocity= Achilles.rawVar.vel(Param.trigAchOn);
    %     EMG(i).isReference = OutputData.TrialSettings.isReference;
    EMGfilt{i}=EMG(i).sync(:,muscleChannels);

    ids = [round(41*fs),round(80*fs)]; % ids that split the signal (in total 18: 9 dynamic FE and 9 RFE)
    torque{i} = EMG(i).torque(ids(1):ids(2));
    angle{i} = EMG(i).pos(ids(1):ids(2));

% time vector

    time{i} = linspace(0,length(EMG(i).sync)/fs,length(EMG(i).sync));
    t = transpose(time{i});

figure, plot(t,EMGfilt{i});
xlabel('time')
ylabel('EMG')

% Remove channel with drastic change

if muscleChannels(1) == 1 || muscleChannels(1) == 65

    [r,c] = size(EMGfilt{i});
     
    for l = 1:c
        for k = 1:r
            if EMGfilt{i}(k,l) == 0
                EMGfilt{i}(:,l) = NaN;
            end
        end
    end
    
    EMGfilt{i} = EMGfilt{i}(:, all(~isnan(EMGfilt{i}), 1));
end

figure, plot(t,EMGfilt{i});
xlabel('time')
ylabel('EMG')

  figure,
  subplot(2,1,1)
    plot(t(ids(1):ids(2)),angle{i});
    hold on
%     xline([40 80],'--k');
    xlabel('time (S)')
    ylabel('Position (rad)')
    subplot(2,1,2)
    plot(t(ids(1):ids(2)),torque{i});
    hold on
%     xline([40 80],'--k');
    xlabel('time (S)')
    ylabel('Torque (Nm)')
end

%% Parallel decomposition

% Checking if a parallel pool is open and if so, shutting it down
% p = gcp('nocreate');
% if isempty(p)        % check if a parallel pool is open
%     parpool(3);
% end

for i=1:length(EMG)  % parfor

    % Notch filter for both decompositions

    for w = 50 %50:50:500  % artifacts at these frequencies
        wo = w/(fs/2);  bw = wo/100;
        [b,a] = iirnotch(wo,bw);
        EMGfilt_notch{i} = filtfilt(b,a,EMGfilt{i});
    end

    [PSD1{i},f1{i}]=period(EMGfilt_notch{i},fs);
    figure, plot(f1{i},PSD1{i});  xlabel('f'); ylabel('PSD'); title('Periodogram');

    % band-pass filter

    fc=[LP,HP]; % useful signal from 10 Hz onwards
    [b,a] = butter(4,fc/(fs/2),'bandpass');
    EMGfilt_band{i} = filtfilt(b,a,EMGfilt_notch{i});

    [PSD2{i},f2{i}]=period(EMGfilt_band{i},fs);

    %     figure, plot(t,EMGfilt_band{i});
    %     legend('40','80')
    [r,c] = size(PSD2{i});

    for l = 1:c-1
        for k = 1:r
            if PSD2{i}(k,l+1)-PSD2{i}(k,l) > 200
                EMGfilt_band{i}(:,l+1) = NaN;
                PSD2{i}(:,l+1) = NaN;
            end
        end
    end

    EMGfilt_band{i} = EMGfilt_band{i}(:, all(~isnan(EMGfilt_band{i}), 1));
    PSD2{i} = PSD2{i}(:, all(~isnan(PSD2{i}), 1));
    figure, plot(f2{i},PSD2{i});  xlabel('f'); ylabel('PSD'); title('Periodogram');

    % Rectification & Envelope

    EMGfilt_rect{i} = abs(EMGfilt_band{i});
    figure, plot(t,EMGfilt_rect{i});
    xlabel('t'); ylabel('Rectified EMG'); title('Rectification');

    fc=5; % useful signal from 10 Hz onwards
    [b,a] = butter(4,fc/(fs/2),'low');
    EMGfilt_Env{i} = filtfilt(b,a,EMGfilt_rect{i});
    figure, plot(t,EMGfilt_Env{i});
    xlabel('t'); ylabel('EMG Envelope'); title('Envelope');

    %         Decomposition

    [ActIndWhole{i},ActIndResid{i},MUPulses{i},Cost{i},MUPulses2{i},Cost2{i},IPTs{i},IPTs2{i},ProcTime{i}] = ...
        GUICKCgrad5(EMGfilt_band{i}(ids(1):ids(2),:)',exFactor,fs,NITER,costFunc,0,iterPlots);      %BSS technique
    tic
    MUPulses{i} = {MUPulses{i}{:} MUPulses2{i}{:}};
    Cost{i} = {Cost{i}{:} Cost2{i}{:}};
    IPTs{i}=[IPTs{i}; IPTs2{i}];
    [MUPulses{i},newInd{i}]=eliminateDuplicateIPT(MUPulses{i},fs,1,0.3);
    IPTs{i} = IPTs{i}(newInd{i},:);
    ProcTime{i} = ProcTime{i} + toc;
end


     %% save/load data after decomposition
    int = 'INTER';
%   save(strcat(path,strcat(trial,num2str(i),mus,int,'.mat')),'MUPulses','IPTs');
    load(strcat(path,strcat(trial,num2str(i),mus,int,'.mat')));

     %% MVC Trial
        
    EMG_norm = cell(1,nRep);
    MVC_trans = EMG_MVC;
    MVC = transpose(MVC_trans);
 
    % Time Vector
    time_MVC = linspace(0,length(MVC)/fs,length(MVC));
    t_MVC = transpose(time_MVC);

    % Notch filter
    w = 50; %50:50:500  % artifacts at these frequencies
    wo = w/(fs/2);  bw = wo/100;
    [b,a] = iirnotch(wo,bw);
    MVC_notch = filtfilt(b,a,MVC);

    [PSD5,f5]=period(MVC_notch,fs);

    % Band Pass Filter
    fc=[LP,HP]; % useful signal from 10 Hz onwards
    [b,a] = butter(4,fc/(fs/2),'bandpass');
    MVCfilt_band = filtfilt(b,a,MVC_notch);

%     [PSD3,f3]=period(MVC,fs);
    [PSD4,f4]=period(MVCfilt_band,fs);

    figure, plot(t_MVC,MVCfilt_band);
    xlabel('t'); ylabel('Bandpass MVC'); title('Filtered MVC');
    figure, plot(t_MVC,MVC);
    xlabel('t'); ylabel('MVC'); title('Actual MVC');

    % Rectification & Envelope
    MVCfilt_rect = abs(MVCfilt_band);
    figure, plot(t_MVC,MVCfilt_rect);
    xlabel('t'); ylabel('Rectified MVC'); title('Rectification');

    fc=5; % useful signal from 10 Hz onwards
    [b,a] = butter(4,fc/(fs/2),'low');
    MVCfilt_Env = filtfilt(b,a,MVCfilt_rect);

%         [PSD5,f5]=period(MVCfilt_Env,fs);

    figure, plot(t_MVC,MVCfilt_Env);
    xlabel('t'); ylabel('MVC Envelope'); title('Envelope');

    %% Mean of MVC trial

    for j = 1:length(MVCfilt_Env)
        mean_MVC(j) = mean(MVCfilt_Env(j,:));
    end
    
        mean_MVC_norm = transpose(mean_MVC);

    [PSD8,f8]=period(mean_MVC_norm,fs);

    figure, plot(t_MVC,mean_MVC_norm);
    xlabel('t'); ylabel('Mean of Normalized MVC'); title('All EMG channels');

    maxi = max(max(mean_MVC_norm));

    %% Mean of all EMG channels

    MU_time_pos = {};
    RT = {};
    RT_final = {};

for i = 1:length(EMGfilt_Env)
    for j = 1:length(EMGfilt_Env{i})
        mean_EMG{i}(j) = mean(EMGfilt_Env{i}(j,:));
    end
        mean_EMG_norm{i} = transpose(mean_EMG{i});

    [PSD7{i},f7{i}]=period(mean_EMG_norm{i},fs);

    figure, plot(t(ids(1):ids(2)),mean_EMG_norm{i}(ids(1):ids(2)));
    hold on
    yyaxis("right")
    plot(t(ids(1):ids(2)),torque{i})
    legend('Mean normalized EMG','Torque')
    xlabel('t'); ylabel('Torque'); title('Mean of all EMG channels');
end

%% Normalizing data - AD of EMG Envelope
   
    for i = 1:nRep
        EMGfilt_mat = mean_EMG_norm{i}(ids(1):ids(2),:)';
        EMG_norm_ind = [];

        for j = 1:length(EMGfilt_mat)
            EMG_norm_ind(:,j) = EMGfilt_mat(:,j)./maxi;
        end
        EMG_norm{i} = EMG_norm_ind';
            
        figure, plot(t(ids(1):ids(2)),EMG_norm{i});
        hold on
    yyaxis("right")
    plot(t(ids(1):ids(2)),torque{i})
    xlabel('t'); ylabel('Mean of all EMG channels'); title('All EMG channels');
%         xlabel('t'); ylabel('Normalized EMG'); title('Normalization');
    
%         figure, plot(t(ids(1):ids(2)),EMGfilt_Env{i}(ids(1):ids(2),:));
%         xlabel('t'); ylabel('EMG Envelope'); title('Envelope without normalization');    
    end

%% save data of bipolar EMG (Bipolar only has EMG Envelope data, so save data here)

% for i = 1:nRep
%     max_EMG_norm_save(i) = max(max(EMG_norm{i}));
%     mean_EMG_norm_save = EMG_norm{i};
%     torque_save = torque{i};
%     angle_save = angle{i};
%     save(strcat(path,strcat(trial,num2str(i),mus,'.mat')),'fs','exFactor','LP','HP','NITER','costFunc','mean_EMG_norm_save','max_EMG_norm_save','torque_save','angle_save');
% end

%% Quality Control - PNR

for i = 1:length(MUPulses)
    for j=1:length(MUPulses{i})
        t_hat{i} = IPTs{i}(j,:)';
        Discharge{i} = MUPulses{i}{j};
        t_hat_1{i} = 1:length(t_hat{i});
        noDischarge{i} = setdiff(t_hat_1{i}',Discharge{i});
        PNR{i}(j) = 10*log10(mean(t_hat{i}(Discharge{i}).^2)/mean(t_hat{i}(noDischarge{i}).^2));
    end
     
    % Visualize PNR

    t_hat_Discharge = cell(1,nRep);
    t_hat_Discharge{i} = NaN(size(t_hat{i}));
    for k = 1:length(t_hat{i})
        for n = 1:length(Discharge{i})
            if k == Discharge{i}(n)
                t_hat_Discharge{i}(k) = t_hat{i}(Discharge{i}(n));
            end
        end
    end

    figure
    plot(t_hat{i})
    hold on
    plot(t_hat_Discharge{i},'o')
    title('PNR with IPT');
end

 %% Omitting spike trains

%  Omitting less number of MU pulses

for i = 1:length(MUPulses)
    for j = 1:length(MUPulses{i})
        if length(MUPulses{i}{j})>500
            MUPulses{i}{j} = [];
            IPTs{i}(j,:) = NaN;
            PNR{i}(j) = NaN;
        elseif length(MUPulses{i}{j})<100
            MUPulses{i}{j} = [];
            IPTs{i}(j,:) = NaN;
            PNR{i}(j) = NaN;
        else
            MUPulses{i}{j} = MUPulses{i}{j};
            IPTs{i}(j,:) = IPTs{i}(j,:);
            PNR{i}(j) = PNR{i}(j);
        end
    end
 
%  Omitting low PNR values

    for j=1:length(MUPulses{i})
        if PNR{i}(j) < 20
            MUPulses{i}{j} = [];
            IPTs{i}(j,:) = NaN;
            PNR{i}(j) = NaN;
        else
            MUPulses{i}{j} = MUPulses{i}{j};
            IPTs{i}(j,:) = IPTs{i}(j,:);
            PNR{i}(j) = PNR{i}(j);
        end
    end
        MUPulses{i} = MUPulses{i}(~cellfun('isempty',MUPulses{i}));
        PNR{i} = PNR{i}(~isnan(PNR{i}));
        IPTs{i} = rmmissing(IPTs{i});
end
        IPTs = IPTs(~cellfun('isempty',(IPTs)));
        MUPulses = MUPulses(~cellfun('isempty',MUPulses));
        PNR = PNR(~cellfun('isempty',PNR));

%% Spike Trains

for i = 1:length(MUPulses)
    plotFlag = 1;
    figure
    [MUDP{i},DR{i},ax] = getSpikeTrains(MUPulses{i}, t(ids(1):ids(2)), fs,torque{i}, plotFlag);
    
    meanDR{i} = mean(DR{i},'omitnan');
    stdDR{i} = std(DR{i},'omitnan');
    IDI{i} = 1./DR{i};
    averageIDI{i} = mean(IDI{i},'omitnan');
    meanstdDR{i} = mean(stdDR{i});
end

%% Omitting spike trains that lie on different flexion

if muscleChannels(1) == 1           %TA -positive cycle (DF)
    for i = 1:length(MUDP)
        for j = 1:length(MUDP{i})
            if torque{i}(j)<0            
                MUDP{i}(j,:) = 0;
                DR{i}(j,:) = NaN;
            else
                MUDP{i}(j,:) = MUDP{i}(j,:);
                DR{i}(j,:) = DR{i}(j,:);
            end
        end
    
        [r2, c2] = size(MUDP{i});
        for j = 1:c2    
            Nz{i}(j) = nnz(MUDP{i}(:,j));
            if MUDP{i}(:,j) == 0
                MUPulses{i}{j} = [];
                IPTs{i}(j,:) = NaN;
                PNR{i}(j) = NaN;
            elseif Nz{i}(j)<25
                MUPulses{i}{j} = [];
                IPTs{i}(j,:) = NaN;
                PNR{i}(j) = NaN;
            end 
            if isempty(MUPulses{i}{j})
                MUDP{i}(:,j) = NaN;
                DR{i}(:,j) = NaN;
                meanDR{i}(j) = NaN;
            end
        end
    
            MUDP{i} = transpose(MUDP{i});
            MUDP{i} = rmmissing(MUDP{i});
            MUDP{i} = transpose(MUDP{i});
            meanDR{i} = meanDR{i}(~isnan(meanDR{i}));
            MUPulses{i} = MUPulses{i}(~cellfun('isempty',MUPulses{i}));
            PNR{i} = PNR{i}(~isnan(PNR{i}));
            IPTs{i} = rmmissing(IPTs{i});
    
    end
            IPTs = IPTs(~cellfun('isempty',(IPTs)));
            MUPulses = MUPulses(~cellfun('isempty',MUPulses));
            PNR = PNR(~cellfun('isempty',PNR));

elseif muscleChannels(1) == 65      %SOL - negative cycle (PF)
    
    for i = 1:length(MUDP)
        for j = 1:length(MUDP{i})
            if torque{i}(j)>0            
                MUDP{i}(j,:) = 0;
                DR{i}(j,:) = NaN;
            else
                MUDP{i}(j,:) = MUDP{i}(j,:);
                DR{i}(j,:) = DR{i}(j,:);
            end
        end
    
        [r2, c2] = size(MUDP{i});
        for j = 1:c2    
            Nz{i}(j) = nnz(MUDP{i}(:,j));
            if MUDP{i}(:,j) == 0
                MUPulses{i}{j} = [];
                IPTs{i}(j,:) = NaN;
                PNR{i}(j) = NaN;
            elseif Nz{i}(j)<25
                MUPulses{i}{j} = [];
                IPTs{i}(j,:) = NaN;
                PNR{i}(j) = NaN;
            end 
            if isempty(MUPulses{i}{j})
                MUDP{i}(:,j) = NaN;
                DR{i}(:,j) = NaN;
                meanDR{i}(j) = NaN;
            end
        end
    
            MUDP{i} = transpose(MUDP{i});
            MUDP{i} = rmmissing(MUDP{i});
            MUDP{i} = transpose(MUDP{i});
            meanDR{i} = meanDR{i}(~isnan(meanDR{i}));
            MUPulses{i} = MUPulses{i}(~cellfun('isempty',MUPulses{i}));
            PNR{i} = PNR{i}(~isnan(PNR{i}));
            IPTs{i} = rmmissing(IPTs{i});
    
    end
            IPTs = IPTs(~cellfun('isempty',(IPTs)));
            MUPulses = MUPulses(~cellfun('isempty',MUPulses));
            PNR = PNR(~cellfun('isempty',PNR));
end

%% Recruitment Threshold

    RT_plot = cell(1,length(MUPulses));

for i = 1:length(MUPulses)
    final_array{i} = first_elements(MUPulses{1,i});
    RT_plot{i} = NaN(size(EMG_norm{i}));
    for j = 1:length(final_array{i})
        for k = 1:length(final_array{i}{1,j})
            [~,MU_time_pos{1,i}{1,j}(k)] = ismember(t(final_array{1,i}{1,j}(k)),t);
            RT{i}{j}(k) = EMG_norm{i}(MU_time_pos{1,i}{1,j}(k));
        end
    end
    
    l=1;

    for j = 1:length(RT{i})
        for k = 1:length(RT{i}{j})
            MU_pos_final{i}(l) = MU_time_pos{i}{j}(k);
            RT_final{i}(l) =  RT{i}{j}(k);
            l=l+1;
        end
    end
    
    for j = 1:length(RT_plot{i})
       for k = 1:length(MU_pos_final{i})
           if j == MU_pos_final{i}(k)
              if RT_final{i}(k)>0
                 RT_plot{i}(j) = RT_final{i}(k);
              else
                 RT_plot{i}(j) = NaN;
              end
           end
       end
    end

   figure,
   plot(RT_final{i},'o');
   title('Recruitment Threshold');

    % Mean RT

    for j = 1:length(final_array{i})
        mean_RT{i}(j) = mean(RT{i}{j});
    end
    figure
    plot(mean_RT{i},'o')
    title('Mean Recruitment Threshold');

    % Normalize DR

    for j = 1:length(MUPulses{i})
        DR_norm{i}(j) = meanDR{i}(j)./40;
    end

    % Data Set for Principal Component Analysis

    Original_dataset{i} = [mean_RT{i}; DR_norm{i}];
    Original_dataset_T{i} = transpose(Original_dataset{i});
end

%% Saving the data of motor units and EMG Envelopes of monopolar EMG

for i = 1:nRep
    max_EMG_norm_save(i) = max(max(EMG_norm{i}));
    MUPulses_save = MUPulses{i};
    EMG_norm_save = EMG_norm{i};
    IPTs_save = IPTs{i};
    MUDP_save = MUDP{i};
    meanDR_save = meanDR{i};
    DR_save = DR{i};
    PNR_save = PNR{i};
    RT_final_save = RT_final{i};
    RT_plot_save = RT_plot{i};
    Original_dataset_T_save = Original_dataset_T{i};
    torque_save = torque{i};
    angle_save = angle{i};
    save(strcat(path,strcat(trial,num2str(i),mus,'.mat')),'fs','exFactor','LP','HP','NITER','costFunc','MUPulses_save','IPTs_save','MUDP_save','DR_save','meanDR_save','PNR_save','RT_final_save','RT_plot_save','Original_dataset_T_save','EMG_norm_save','max_EMG_norm_save','torque_save','angle_save');
end
