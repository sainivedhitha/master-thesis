clc
clear all
close all

%% Import data

directory = 'D:\SaiArunajatesan\Thesis\Multisegment algorithm - Stiffness via system identification\Simulation files\Simulation Results\';
data = {'trial_1-','trial_2-','trial_3-','trial_4-','trial_5-','trial_6-'};
sub = {'subject01\','subject02\','subject03\','subject04\','subject05\','subject06\'};
static = 'Static-';
path1 = 'D:\SaiArunajatesan\Thesis\Data\';
path = 'D:\SaiArunajatesan\Thesis\Data\subject01\';

muscles = {'TA';'GM';'GL';'PL'};
nRep = 8;
fs = 2048;
ids = [round(41*fs),round(80*fs)];

% time = linspace(41,80,length(ids(1):ids(2)));
time = linspace(45,47,1000);

t = transpose(time);

%% load stiffness results

trial = 'Trial5C1R';

for i = 1:6
    Output_static(i) = load(strcat(path1, strcat(sub{i}, strcat(static, num2str(i)))));        
    stiff_exp{i} = Output_static(i).stiffness;
    for j = 1:8
        Output_simulation{i}(j) = load(strcat(directory, strcat(data{i}, num2str(j))));
        stiff_sim{i}{j} = Output_simulation{i}(j).stiffness;
    end
end

Output_simulation{4}(9) = load(strcat(directory, strcat(data{4}, num2str(9))));
stiff_sim{4}{9} = Output_simulation{4}(9).stiffness;

%% load Activation Dynamics

Output_AD = load(strcat(path,'PCA_TA','.mat'));

AD = Output_AD.AD;
    
%% load Normalized EMG

 for i = 1:nRep   
    for j=1:length(muscles)
        if j == 1
            EMG_norm{i}(j) = load(strcat(path,strcat(trial,num2str(i),muscles{j},'.mat')));
        else
            EMG_norm1{i}(j) = load(strcat(path,strcat(trial,num2str(i),muscles{j},'.mat')));
        end
    end
 end
 
%% Moving median value of stiffness through SI

for j = 1:6
    stiff_exp{1,j} = movmedian(stiff_exp{1,j},5);
end

for j = 1:6
    stiff_exp{1,j} = movmean(stiff_exp{1,j},5);
end
%% Save filtered stiffness data

for i =1:6
    stiffness = stiff_exp{i};
    save(strcat('D:\SaiArunajatesan\Thesis\Data\',num2str(i),'.mat'),'stiffness');
end

%% Plot angle,torque & stiffness

for i = 1:nRep
    angle{i} = 10+(EMG_norm1{i}(2).angle_save-0.9447)*pi/180;
    torque{i} = EMG_norm1{i}(2).torque_save;
end

for i = 1:nRep
       figure,
        ax1 = subplot(311);
        plot(angle{i})
%         ylim([9 10])
        title('Angle(deg)');
        hold on
        ax2 = subplot(312);
        plot(torque{i})
%         ylim([15 40])
        title('Torque');
        hold on
        ax3 = subplot(313);
        plot(stiff_exp{3}(:,10))
%         ylim([15 40])
        title('Stiffness');
        hold off
        linkaxes([ax1 ax2 ax3],'x')
end

%% Segmentation of torque & angle

dt = 1/fs;                                          % Sampling time [s]
f = 0.5;                                            % Sine frequency [Hz]
minPeakHeight = 2.5;

for i = 1:nRep
    [pks{i},locs{i}] = findpeaks(-torque{i}, 'MinPeakHeight', minPeakHeight,'MinPeakDistance', 0.7*ceil(fs/f));
    segments{i} = diff(locs{i});
end

for i = 1:length(locs)
    for j = 1
        torque_seg{i}{j} = torque{i}(1:locs{i}(j),1);
        angle_seg{i}{j} = angle{i}(1:locs{i}(j),1);
    end
end

for i = 1:length(locs)
    for j = 2:length(locs{i})
        torque_seg{i}{j} = torque{i}(locs{i}(j-1)+1:locs{i}(j),1);
        angle_seg{i}{j} = angle{i}(locs{i}(j-1)+1:locs{i}(j),1);
    end
end

for i = 1:length(locs)
    torque_seg{i}{1} = [];
    angle_seg{i}{1} = [];
    torque_seg{i} = torque_seg{i}(~cellfun('isempty',torque_seg{i}));
    angle_seg{i} = angle_seg{i}(~cellfun('isempty',angle_seg{i}));
end

%

for i = 1:length(torque)
    locs_plot{i} = NaN(size(torque{i}));
    for j = 1:length(torque{i})    
        for k = 1:length(pks{i})
            if j == locs{i}(k)
                locs_plot{i}(j) = pks{i}(k);
            end
        end
    end
end


%% Segmentation of EMG envelope

for i = 1:nRep
    EMG_norm_save_TA{i} = EMG_norm{i}.EMG_norm_save;
    EMG_norm_save_GM{i} = EMG_norm1{i}(2).EMG_norm_save;
    EMG_norm_save_GL{i} = EMG_norm1{i}(3).EMG_norm_save;
    EMG_norm_save_PL{i} = EMG_norm1{i}(4).EMG_norm_save;
 end
 
 for j=1:nRep
     for i = 1
        EMG_TA_seg{j}{i} = EMG_norm_save_TA{j}(1:locs{j}(i),1);
        EMG_GM_seg{j}{i} = EMG_norm_save_GM{j}(1:locs{j}(i),1);
        EMG_GL_seg{j}{i} = EMG_norm_save_GL{j}(1:locs{j}(i),1);
        EMG_PL_seg{j}{i} = EMG_norm_save_PL{j}(1:locs{j}(i),1);
    end

    for i = 2:length(locs{j})
        EMG_TA_seg{j}{i} = EMG_norm_save_TA{j}(locs{j}(i-1)+1:locs{j}(i),1);
        EMG_GM_seg{j}{i} = EMG_norm_save_GM{j}(locs{j}(i-1)+1:locs{j}(i),1);
        EMG_GL_seg{j}{i} = EMG_norm_save_GL{j}(locs{j}(i-1)+1:locs{j}(i),1);
        EMG_PL_seg{j}{i} = EMG_norm_save_PL{j}(locs{j}(i-1)+1:locs{j}(i),1);    
    end
 end
 
 for i = 1:length(locs)
    EMG_TA_seg{i}{1} = [];
    EMG_GM_seg{i}{1} = [];
    EMG_TA_seg{i} = EMG_TA_seg{i}(~cellfun('isempty',EMG_TA_seg{i}));
    EMG_GM_seg{i} = EMG_GM_seg{i}(~cellfun('isempty',EMG_GM_seg{i}));
    EMG_GL_seg{i}{1} = [];
    EMG_PL_seg{i}{1} = [];
    EMG_GL_seg{i} = EMG_GL_seg{i}(~cellfun('isempty',EMG_GL_seg{i}));
    EMG_PL_seg{i} = EMG_PL_seg{i}(~cellfun('isempty',EMG_PL_seg{i}));
 end
 
%% Segmentation of AD

 for j=1:nRep
    for i = 1
        AD_seg{j}{i} = AD{j}(1:locs{j}(i),1);
    end

    for i = 2:length(locs{j})
        AD_seg{j}{i} = AD{j}(locs{j}(i-1)+1:locs{j}(i),1);
    end
 end
 
 for i = 1:length(locs)
    AD_seg{i}{1} = [];
    AD_seg{i} = AD_seg{i}(~cellfun('isempty',AD_seg{i}));
 end

  %% plot segments for checking

 % Torque
 for i =1:nRep
    figure,
    ax1 = subplot(211)
    plot(torque{i});
    hold on
    xline(locs{i},'--k');
    yticks([-8 0 8])
    ax2 = subplot(212)
    plot(JS{i});
    hold on
    xline(locs{i},'--k');
 end

% EMG envelopes
for i =1:nRep
    figure,
    ax1 = subplot(211)
    plot(EMG_norm_save_GL{i});
    hold on
    xline(locs{i},'--k');
%     yticks([-8 0 8])
    ax2 = subplot(212)
    plot(torque{i});
    hold on
    xline(locs{i},'--k');
end

for i = 1:nRep
    figure,
        plot(EMG_norm{i}.EMG_norm_save);
        hold on
        plot(EMG_norm1{i}(2).EMG_norm_save);
        hold on
        plot(EMG_norm1{i}(3).EMG_norm_save);
        hold on
        plot(EMG_norm1{i}(4).EMG_norm_save);
        hold on
        yyaxis("right")
        plot(torque{i});
      legend('TA','GM','GL','PL','torque')
end

% AD of motor units

for i = 1:nRep
    figure,
        plot(AD{i});
        hold on
        yyaxis("right")
        plot(torque{i});
      legend('AD','torque')
end

for i =1:nRep
    figure,
    ax1 = subplot(211)
    plot(AD{i});
    hold on
    xline(locs{i},'--k');
%     yticks([-8 0 8])
    ax2 = subplot(212)
    plot(torque{i});
    hold on
    xline(locs{i},'--k');
end

%% Time Normalization

% Torque & Angle

 for i = 1:length(torque_seg)
    for j = 1:length(torque_seg{i})
        idx = length(torque_seg{i}{j});
        tor_pad{i}{j} = padarray(torque_seg{i}{j},round(idx*0.3),'replicate','both');
        ang_pad{i}{j} = padarray(angle_seg{i}{j},round(idx*0.3),'replicate','both');
        tor_1600{i}{j}  = resample(tor_pad{i}{j},1600,length(tor_pad{i}{j}));
        ang_1600{i}{j}  = resample(ang_pad{i}{j},1600,length(ang_pad{i}{j}));
        tor_1000{i}{j} = tor_1600{i}{j}(301:1300);
        ang_1000{i}{j} = ang_1600{i}{j}(301:1300);
    end
 end

% Stiffness

    idx = length(stiff_exp{3}(:,10));                                                                    % Change for each subject
    stiff_pad = padarray(stiff_exp{3}(:,10),round(idx*0.3),'replicate','both');                          % Change for each subject
    stiff_1600  = resample(stiff_pad,1600,length(stiff_pad));
    stiff_1000 = stiff_1600(301:1300);

% EMG Envelope

 for i = 1:nRep
    for j = 1:length(EMG_TA_seg{i})
        idx = length(EMG_TA_seg{i}{j});
        EMG_TA_pad{i}{j} = padarray(EMG_TA_seg{i}{j},round(idx*0.3),'replicate','both');
        EMG_GM_pad{i}{j} = padarray(EMG_GM_seg{i}{j},round(idx*0.3),'replicate','both');
        EMG_GL_pad{i}{j} = padarray(EMG_GL_seg{i}{j},round(idx*0.3),'replicate','both');
        EMG_PL_pad{i}{j} = padarray(EMG_PL_seg{i}{j},round(idx*0.3),'replicate','both');
        EMG_TA_1600{i}{j}  = resample(EMG_TA_pad{i}{j},1600,length(EMG_TA_pad{i}{j}));
        EMG_GM_1600{i}{j}  = resample(EMG_GM_pad{i}{j},1600,length(EMG_GM_pad{i}{j}));
        EMG_GL_1600{i}{j}  = resample(EMG_GL_pad{i}{j},1600,length(EMG_GL_pad{i}{j}));
        EMG_PL_1600{i}{j}  = resample(EMG_PL_pad{i}{j},1600,length(EMG_PL_pad{i}{j}));
        EMG_TA_1000{i}{j} = EMG_TA_1600{i}{j}(301:1300);
        EMG_GM_1000{i}{j} = EMG_GM_1600{i}{j}(301:1300);
        EMG_GL_1000{i}{j} = EMG_GL_1600{i}{j}(301:1300);
        EMG_PL_1000{i}{j} = EMG_PL_1600{i}{j}(301:1300);
    end
 end

 % Acivation Dynamics

 for i = 1:nRep
    for j = 1:length(AD_seg{i})
        idx = length(AD_seg{i}{j});
        AD_pad{i}{j} = padarray(AD_seg{i}{j},round(idx*0.3),'replicate','both');
        AD_1600{i}{j}  = resample(AD_pad{i}{j},1600,length(AD_pad{i}{j}));
        AD_1000{i}{j} = AD_1600{i}{j}(301:1300);
    end
 end

 %% Plot torque data

for i = 1:length(EMG_TA_seg{1})
    figure,
    plot(linspace(0, 1000, length(EMG_TA_seg{1}{i})), EMG_TA_seg{1}{i})
    hold on
    plot(linspace(0,1000, length(EMG_TA_pad{1}{i})), EMG_TA_pad{1}{i})
    hold on
    plot(linspace(0, 1000, length(EMG_TA_1600{1}{i})), EMG_TA_1600{1}{i})
    hold on
    plot(linspace(0,1000, length(EMG_TA_1000{1}{i})), EMG_TA_1000{1}{i})
    legend('originalSignal','paddedSignal','resampledSignal','resampled-afterRemovingpadding')
end

% Stiffness data

for i = 1
    figure,
    plot(linspace(0,1600,length(stiff_exp{3}(:,10))),stiff_exp{3}(:,10))
    hold on
    plot(linspace(0,1600,length(stiff_1600)),stiff_1600)
    legend('originalSignal','resampledSignal')
end

 %% Mean of all segments
 
 % Torque & Angle

for i = 1:nRep
    tor = cell2mat(tor_1000{i});
    ang = cell2mat(ang_1000{i});
    tor_1000_sub{i} = tor;
    mean_tor_1000{i} = mean(tor_1000_sub{i},2);
    ang_1000_sub{i} = ang;
    mean_ang_1000{i} = mean(ang_1000_sub{i},2);
end

torq_final_seg = cell2mat(mean_tor_1000);
torq_final = mean(torq_final_seg,2);
angle_final_seg = cell2mat(mean_ang_1000);
angle_final = mean(angle_final_seg,2);

% EMG Envelope

for i = 1:nRep
    EMG_TA = cell2mat(EMG_TA_1000{i});
    EMG_GM = cell2mat(EMG_GM_1000{i});
    TA_1000_sub{i} = EMG_TA;
    mean_TA_1000{i} = mean(TA_1000_sub{i},2);
    GM_1000_sub{i} = EMG_GM;
    mean_GM_1000{i} = mean(GM_1000_sub{i},2);
    EMG_GL = cell2mat(EMG_GL_1000{i});
    EMG_PL = cell2mat(EMG_PL_1000{i});
    GL_1000_sub{i} = EMG_GL;
    mean_GL_1000{i} = mean(GL_1000_sub{i},2);
    PL_1000_sub{i} = EMG_PL;
    mean_PL_1000{i} = mean(PL_1000_sub{i},2);
end

TA_final_seg = cell2mat(mean_TA_1000);
EMG_norm_TA = mean(TA_final_seg,2);
GM_final_seg = cell2mat(mean_GM_1000);
EMG_norm_GM = mean(GM_final_seg,2);
GL_final_seg = cell2mat(mean_GL_1000);
EMG_norm_GL = mean(GL_final_seg,2);
PL_final_seg = cell2mat(mean_PL_1000);
EMG_norm_PL = mean(PL_final_seg,2);

% AD

for i = 1:nRep
    AD_TA = cell2mat(AD_1000{i});
    AD_TA_sub{i} = AD_TA;
    mean_AD_TA{i} = mean(AD_TA_sub{i},2);
end

AD_TA_final_seg = cell2mat(mean_AD_TA);
AD_mean_final = mean(AD_TA_final_seg,2);

%% 
figure,
    plot(t,EMG_norm_TA);
    hold on
    plot(t,EMG_norm_GM);
    hold on
    plot(t,EMG_norm_GL);
    hold on
    plot(t,EMG_norm_PL);
    yyaxis("left")
    ylabel('Activation Dynamics')
    hold on
    yyaxis("right")
    plot(t,torq_final);
    xlabel('time (seconds)'); 
    ylabel('Torque (Nm)');
    legend('TA','GM','GL','PL')
    title('AD of EMG Envelope ')
        
 figure,
    plot(t,AD_mean_final);
    yyaxis("left")
    ylabel('Activation Dynamics')
    hold on
    yyaxis("right")
    plot(t,torq_final);
    ylabel('Torque (Nm)');
    xlabel('time (seconds)'); 

    
    %% Save results for TA & Soleus 
    
    % Here, we omit soleus for motor units, so only it's EMG Envelopes were
    % saved to check the results. But they were worse at the end.
    
%      EMG_norm_SOL = EMG_norm_TA;
%     save(strcat(path,'Results_SOL.mat'),'EMG_norm_SOL');
    save(strcat(path,'Results.mat'),'AD_mean_final','torq_final','EMG_norm_PL','EMG_norm_GL','EMG_norm_GM','EMG_norm_TA','stiff_exp','stiff_1000');
     
%% Save data for calibration (2 seconds / length(1000) data)
        
 % Save torque for calibration

    L = length(torq_final);
    t1 = t;
    knee_torque = zeros(L,1);
    data = [t1,knee_torque,torq_final];
    [r,c] = size(data);
    temp_fn = ['ID' '.sto'];
    fileID = fopen(temp_fn,'w');
    fprintf(fileID,'%s\n','Inverse Dynamics Generalized Forces');
    fprintf(fileID,'%s','Version=');
    fprintf(fileID,'%1.0f\n',1);
    fprintf(fileID,'%s','nRows=');
    fprintf(fileID,'%4.0f\n',r);
    fprintf(fileID,'%s','nColumns=');
    fprintf(fileID,'%1.0f',c);
    fprintf(fileID,'\n%s','inDegrees=no');
    fprintf(fileID,'\n%s','endheader');
    fprintf(fileID,'\n%4s  \t%12s\t%12s\t\n','time','knee_angle_r_moment','ankle_angle_r_moment');
    fprintf(fileID,'%4.4f\t%4.6f\t%4.6f\t\n',data');
    fclose(fileID);

    % Save angle for calibration

    L = length(angle_final);
    t1 = t;
    knee_angle = zeros(L,1)-30;
    data = [t1,knee_angle,angle_final];
    [r,c] = size(data);
    temp_fn = ['IK' '.sto'];
    fileID = fopen(temp_fn,'w');
    fprintf(fileID,'%s\n','Results');
    fprintf(fileID,'%s','Version=');
    fprintf(fileID,'%1.0f\n',1);
    fprintf(fileID,'%s','nRows=');
    fprintf(fileID,'%4.0f\n',r);
    fprintf(fileID,'%s','nColumns=');
    fprintf(fileID,'%1.0f',c);
    fprintf(fileID,'\n%s','inDegrees=yes');
    fprintf(fileID,'\n%s','endheader');
    fprintf(fileID,'\n%4s  \t%12s\t%12s\t\n','time','knee_angle_r','ankle_angle_r');
    fprintf(fileID,'%4.4f\t%4.6f\t%4.6f\t\n',data');
    fclose(fileID);

    % Save stiffness for calibration

    L = length(stiff_1000);
    t1 = t;
    knee_stiff = zeros(L,1);
    data = [t1,knee_stiff,stiff_1000];
        [r,c] = size(data);
    temp_fn = ['jointStiffness' '.sto'];
    fileID = fopen(temp_fn,'w');
    fprintf(fileID,'%s\n','Inverse Dynamics Generalized Forces');
    fprintf(fileID,'%s','Version=');
    fprintf(fileID,'%1.0f\n',1);
    fprintf(fileID,'%s','nRows=');
    fprintf(fileID,'%4.0f\n',r);
    fprintf(fileID,'%s','nColumns=');
    fprintf(fileID,'%1.0f',c);
    fprintf(fileID,'\n%s','inDegrees=no');
    fprintf(fileID,'\n%s','endheader');
    fprintf(fileID,'\n%4s  \t%12s\t%12s\t\n','time','knee_angle_r_moment','ankle_angle_r_moment');
    fprintf(fileID,'%4.4f\t%4.6f\t%4.6f\t\n',data');
    
    % Save AD for calibration

    t1 = t;
    data = [t1,AD_mean_final,EMG_norm_GM,EMG_norm_GL,EMG_norm_PL];        
    [r,c] = size(data);
    temp_fn = ['AD' '.sto'];
    fileID = fopen(temp_fn,'w');
    fprintf(fileID,'%s\n','Activation Dynamics');
    fprintf(fileID,'%s','nRows=');
    fprintf(fileID,'%4.0f\n',r);
    fprintf(fileID,'%s','nColumns=');
    fprintf(fileID,'%1.0f',c);
    fprintf(fileID,'\n%s','endheader');
    fprintf(fileID,'\n%4s  \t%9s\t%8s\t%8s\t%9s\t\n','time','tibant_r','gasmed_r','gaslat_r','perlong_r');
    fprintf(fileID,'%4.4f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t\n',data');
    
     % Save Normalized EMG envelopes for calibration

    t1 = t;
    data = [t1,EMG_norm_TA,EMG_norm_GM,EMG_norm_GL,EMG_norm_PL];        
    [r,c] = size(data);
    temp_fn = ['EMG' '.sto'];
    fileID = fopen(temp_fn,'w');
    fprintf(fileID,'%s\n','Normalized EMG Linear Envelopes');
    fprintf(fileID,'%s','nRows=');
    fprintf(fileID,'%4.0f\n',r);
    fprintf(fileID,'%s','nColumns=');
    fprintf(fileID,'%1.0f',c);
    fprintf(fileID,'\n%s','endheader');
    fprintf(fileID,'\n%4s  \t%9s\t%8s\t%8s\t%9s\t\n','time','tibant_r','gasmed_r','gaslat_r','perlong_r');
    fprintf(fileID,'%4.4f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t\n',data');
 
    %% Save Data for execution (entire 40 seconds data)
    
    time = linspace(0,245790/fs,245790);
    t = transpose(time);
    
    % Save normalized EMG envelopes
    
for i = 1:nRep    
    L = length(ids(1):ids(2));
    t1 = t(ids(1):ids(2));
    data = [t1,EMG_norm{i}.EMG_norm_save,EMG_norm1{i}(2).EMG_norm_save,EMG_norm1{i}(3).EMG_norm_save,EMG_norm1{i}(4).EMG_norm_save];
    [r,c] = size(data);
    temp_fn = ['emgFilt' '_' num2str(i) '.sto'];
    fileID = fopen(temp_fn,'w');
    fprintf(fileID,'%s\n','Normalized EMG Linear Envelopes');
    fprintf(fileID,'%s','nRows=');
    fprintf(fileID,'%4.0f\n',r);
    fprintf(fileID,'%s','nColumns=');
    fprintf(fileID,'%1.0f',c);
    fprintf(fileID,'\n%s','endheader');
    fprintf(fileID,'\n%4s  \t%9s\t%8s\t%8s\t%9s\t\n','time','tibant_r','gasmed_r','gaslat_r','perlong_r');
    fprintf(fileID,'%4.4f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t\n',data');
    fclose(fileID);
end

    % Save torque
    
for i=1:nRep
    L = length(ids(1):ids(2));
    t1 = t(ids(1):ids(2));
    knee_torque = zeros(L,1);
    data = [t1,knee_torque,torque{i}];
    [r,c] = size(data);
    temp_fn = ['ID' '_' num2str(i) '.sto'];
    fileID = fopen(temp_fn,'w');
    fprintf(fileID,'%s\n','Inverse Dynamics Generalized Forces');
    fprintf(fileID,'%s','Version=');
    fprintf(fileID,'%1.0f\n',1);
    fprintf(fileID,'%s','nRows=');
    fprintf(fileID,'%4.0f\n',r);
    fprintf(fileID,'%s','nColumns=');
    fprintf(fileID,'%1.0f',c);
    fprintf(fileID,'\n%s','inDegrees=no');
    fprintf(fileID,'\n%s','endheader');
    fprintf(fileID,'\n%4s  \t%12s\t%12s\t\n','time','knee_angle_r_moment','ankle_angle_r_moment');
    fprintf(fileID,'%4.4f\t%4.6f\t%4.6f\t\n',data');
    fclose(fileID);
end

    % Save angle
    
for i=1:nRep
    L = length(ids(1):ids(2));
    t1 = t(ids(1):ids(2));
    knee_angle = zeros(L,1)-30;
    data = [t1,knee_angle,angle{i}];
    [r,c] = size(data);
    temp_fn = ['IK' '_' num2str(i) '.sto'];
    fileID = fopen(temp_fn,'w');
    fprintf(fileID,'%s\n','Results');
    fprintf(fileID,'%s','Version=');
    fprintf(fileID,'%1.0f\n',1);
    fprintf(fileID,'%s','nRows=');
    fprintf(fileID,'%4.0f\n',r);
    fprintf(fileID,'%s','nColumns=');
    fprintf(fileID,'%1.0f',c);
    fprintf(fileID,'\n%s','inDegrees=yes');
    fprintf(fileID,'\n%s','endheader');
    fprintf(fileID,'\n%4s  \t%12s\t%12s\t\n','time','knee_angle_r','ankle_angle_r');
    fprintf(fileID,'%4.4f\t%4.6f\t%4.6f\t\n',data');
end

    % Save AD

for i = 1:nRep
    t1 = t(ids(1):ids(2));
    L = length(t1);
    data = [t1,AD{i},EMG_norm1{i}(2).EMG_norm_save,EMG_norm1{i}(3).EMG_norm_save,EMG_norm1{i}(4).EMG_norm_save];
    [r,c] = size(data);
    temp_fn = ['AD' '_' num2str(i) '.sto'];
    fileID = fopen(temp_fn,'w');
    fprintf(fileID,'%s\n','Activation Dynamics');
    fprintf(fileID,'%s','nRows=');
    fprintf(fileID,'%4.0f\n',r);
    fprintf(fileID,'%s','nColumns=');
    fprintf(fileID,'%1.0f',c);
    fprintf(fileID,'\n%s','endheader');
    fprintf(fileID,'\n%4s  \t%9s\t%8s\t%8s\t%9s\t\n','time','tibant_r','gasmed_r','gaslat_r','perlong_r');
    fprintf(fileID,'%4.4f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t\n',data');
end
