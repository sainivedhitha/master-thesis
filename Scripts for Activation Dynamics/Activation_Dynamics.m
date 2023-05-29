% Run until data extraction section in
% parallelDecompositionCKC_MUinhibition

%% Load results of parallel decomposition

for i = 1:nRep
    Output(i) = load(strcat(path,strcat(trial,num2str(i),mus,'.mat')));
end

%% Plot PNR for reference

PNR1 = [];
k=1;

for i = 1:nRep
    PNR{i} = Output(i).PNR_save;    
    PNR_mean(i) = mean(PNR{i}); 
    for j = 1:length(PNR{1,i})
        PNR1(k) =  PNR{i}(j);
        k=k+1;
    end
end

figure
for i =1:nRep
    plot(PNR{i},'o');
    hold on
end
xlabel('No. of MU pulses'); ylabel('PNR');

% Histogram
   
figure,
    histogram(PNR1);
    xlabel('PNR (dB)'); ylabel('Number of MU Pulses')
    ylim([0 20])

%% PCA Analysis

for i = 1:nRep
    PCA_dataset{i} = Output(i).Original_dataset_T_save;
    DischargeRate{i} = Output(i).meanDR_save;
    RecruitmentThreshold{i} = Output(i).RT_final_save;
end

PCA_dataset_final = vertcat(PCA_dataset{:});
DischargeRate_final = cell2mat(DischargeRate);
RecruitmentThreshold_final = cell2mat(RecruitmentThreshold);

% Subtract mean from data
mean_DR_sub = mean(PCA_dataset_final(:,2),'omitnan');
mean_DR_zero = PCA_dataset_final(:,2) - mean_DR_sub;
mean_RT_sub = mean(PCA_dataset_final(:,1));
mean_RT_zero = PCA_dataset_final(:,1) - mean_RT_sub;

figure, 
scatter(PCA_dataset_final(:,1),PCA_dataset_final(:,2),25,'k','filled')
xlabel('RT'); ylabel('DR');

figure, 
scatter(mean_RT_zero,mean_DR_zero,25,'k','filled')
xlabel('RT'); ylabel('DR');
title('RT (zero) versus DR (zero)');

Dataset_zero = [mean_RT_zero, mean_DR_zero];

[PCA, score, latent] = pca(Dataset_zero);

figure,
biplot(PCA(:,1:2),'scores',score(:,1:2),'varlabels',{'v_1','v_2'});
xlabel('PC1(RT)'); ylabel('PC2(Normalized DR)');
title('RT versus norm DR with eigen vector');

if latent(1) >latent(2)
    PC1 = score(:,1);
    PC2 = score(:,2);
else
    PC1 = score(:,2);
    PC2 = score(:,1);
end

var_PC1 = var(PC1);
var_PC2 = var(PC2);

figure, 
plot(PC1, PC2,'ko','MarkerFaceColor','g')
title('Principal Components (PC1 Vs PC2)');

figure, 
biplot(PCA(:,1:2));
hold on
scatter(mean_RT_zero,mean_DR_zero,25,'k','filled')
xlabel('RT'); ylabel('Normalized DR');
title('RT versus norm DR with eigen vector');

figure, 
biplot(PCA(:,1:2));
hold on
scatter(PCA_dataset_final(:,1),PCA_dataset_final(:,2),25,'k','filled')
xlabel('RT'); ylabel('Normalized DR');
title('RT versus norm DR with eigen vector');

%% Linear Mapping 

% size mapping

lim_1 = min(PC1);
lim_2 = max(PC1);
lim_3 = 2+min(PC1);

Trans_map = ((lim_1+lim_3)+(lim_3-lim_1)*((2*PC1 - (lim_1+lim_2))/(lim_2-lim_1)))/2;

figure, 
plot(PC1,'ko','MarkerFaceColor','g')
title('First Principal Component (PC1)');

figure, 
plot(Trans_map,'ko','MarkerFaceColor','b')
title('First Principal Component (PC1) - size mapping');

% Contraction times mapping

lim_4 = min(Trans_map);
lim_5 = max(Trans_map);

if muscleChannels(1) == 1
    lim_7 = 0.134;            % Tc of TA
    lim_6 = 0.047;            % Tc of TA
elseif muscleChannels(1) == 65
    lim_7 = 0.1859;           % Tc of SOL
    lim_6 = 0.1271;           % Tc of SOL
end

Tc_map = ((lim_6+lim_7)+(lim_7-lim_6)*((2*PC1 - (lim_4+lim_5))/(lim_5-lim_4)))/2;

figure, 
plot(Tc_map,'ko','MarkerFaceColor','y')
title('First Principal Component (PC1) - Tc mapping');

% Peak amplitude mapping

lim_8 = 0.1;            
lim_9 = 1;             

Ap_map = ((lim_8+lim_9)+(lim_9-lim_8)*((2*PC1 - (lim_4+lim_5))/(lim_5-lim_4)))/2;

figure, 
plot(Ap_map,'ko','MarkerFaceColor','m')
title('First Principal Component (PC1) - Ap mapping');

%% Apply all the values in the formula (Convolution)

for  i = 1:nRep
    length1(i) = length(Output(i).MUPulses_save);
end

for i = 1
    for j = 1:length1(i)
        Ap_map1{i}(j) = Ap_map(j);
        Tc_map1{i}(j) = Tc_map(j);
    end
end

for i = 2:length(length1)
    for j = 1:length1(i-1)+length1(i)
        if j <= length1(i-1)
            Ap_map1{i}(j) = NaN;
            Tc_map1{i}(j) = NaN;
        else
            Ap_map1{i}(j) = Ap_map(j);
            Tc_map1{i}(j) = Tc_map(j);     
        end
    end
end

for i = 2:length(length1)
            Ap_map1{i} = Ap_map1{i}(~isnan(Ap_map1{i}));
            Tc_map1{i} = Tc_map1{i}(~isnan(Tc_map1{i}));
end

for i = 1:length(Ap_map1)
    for j = 1:length(Ap_map1{i})
        Ap_norm{i}(j) = Ap_map1{i}(j)./length(Ap_map1{i});
    end
end

fs = 2048;
T = 1/fs;
for  i = 1:nRep
    MUDP_data{i} = Output(i).MUDP_save;
end

for i = 1:length(MUDP_data)
    [r,c] = size(MUDP_data{i});
        for j = 1:c
            for n = 1:r
                if n == 1 || n == 2
                    f_n{i}{j}(n) = 0;
                else
                    f_n{i}{j}(n) = 2*exp(-T/Tc_map1{i}(j))*f_n{i}{j}(n-1) - (exp(-2*T/Tc_map1{i}(j))*f_n{i}{j}(n-2))+(Ap_norm{i}(j)*(T^2)/Tc_map1{i}(j))*exp(1-T/Tc_map1{i}(j))*MUDP_data{i}(n-1,j);                 end
            end
        end
end

for i = 1:nRep
    f_n_A.f_n{i} = f_n{i};
    if length(f_n_A.f_n{i}) > 1
        f_n_A.fn{i} = sum(cell2mat(f_n_A.f_n{i}'));
    else
        f_n_A.fn{i} = cell2mat(f_n_A.f_n{i});
    end
    fn{i} = f_n_A.fn{i};
end

%%
max_ii = [];

figure,
for i = 1:nRep
    max_ii = [max_ii,max(fn{i})];
    plot(t(ids(1):ids(2)),fn{i});
    hold on
end
legend('Trial 1','Trial 2','Trial 3','Trial 4','Trial 5','Trial 6','Trial 7','Trial 8');
xlabel('time');
ylabel('f(n)');

figure,
plot(max_ii,'o');

% T - inverse of fs
% f(n-1) - for first two are 0
% e(n) - MUDP
% plot f(n) vs time
% add all f(n) per trial

%% Normalizing f(n)  - AD of motor units

for i = 1:nRep
    max_EMG_norm(i) = max(Output(i).max_EMG_norm_save);
    lim_11 = max(max_EMG_norm);
    EMG_norm{i} = transpose(Output(i).EMG_norm_save);
end

for i = 1:nRep
    lim_10 = min(fn{i});
    lim_12 = min(fn{i});
    lim_13 = max(fn{i});
    AD{i} = transpose(((lim_10+lim_11)+(lim_11-lim_10)*((2*fn{i} - (lim_12+lim_13))/(lim_13-lim_12)))/2);
end

%% plot AD of EMG, Motor units & torque

for i = 1:nRep
    figure
    plot(t(ids(1):ids(2)),AD{i});
    hold on    
    plot(t(ids(1):ids(2)),EMG_norm{i}');
    hold on
    yyaxis("right")
    plot(t(ids(1):ids(2)),torque{i});
    legend('Activation Dynamics','Envelope','torque');
    xlabel('time');
    ylabel('Torque');
end

%% Saving the results

dir=strcat(strcat(path,strcat('PCA_',mus,'.mat')));
save(dir,'PCA_dataset_final','PNR','Dataset_zero','PCA','PC1','PC2','Trans_map','Tc_map','Ap_map','AD')
