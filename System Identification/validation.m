    % Run Stiffness code until the segmentation of torque & angle

%% load data of trials from CEINMS and reference

% ceinmsTestData3 - AD of motor units
% ceinmsTestData1 - AD of EMG Envelope

    path = 'C:\Users\s2662841\Documents\ceinms_rt\cfg\Subject 1\ceinmsTestData3\CalibrationFiles\';
    path_calib = 'C:\Users\s2662841\Documents\ceinms_rt\cfg\Subject 1\ceinmsTestData3\CalibrationFiles\Calibration\';
    path1 = 'C:\Users\s2662841\Documents\ceinms_rt\cfg\Subject 1\ceinmsTestData1\CalibrationFiles\';
    path1_calib = 'C:\Users\s2662841\Documents\ceinms_rt\cfg\Subject 1\ceinmsTestData3\CalibrationFiles\Calibration\';
    result = 'Results\';
    nRep = 8;
    trial = {'Trial1\','Trial2\','Trial3\','Trial4\','Trial5\','Trial6\','Trial7\','Trial8\'};
    load('D:\SaiArunajatesan\Thesis\Data\subject01\Static-1.mat');                                          

    for i = 1:nRep
        ID_AD(i) = importdata(strcat(path,strcat(trial{i},strcat(result,'Torque.sto'))));
        ID_EMG(i) = importdata(strcat(path1,strcat(trial{i},strcat(result,'Torque.sto'))));
        JS_AD(i) = importdata(strcat(path,strcat(trial{i},strcat(result,'JointStiffness.sto'))));
        JS_EMG(i) = importdata(strcat(path1,strcat(trial{i},strcat(result,'JointStiffness.sto'))));
        emgFilt_AD(i) = importdata(strcat(path,strcat(trial{i},strcat(result,'emg.sto'))));
        emgFilt_EMG(i) = importdata(strcat(path1,strcat(trial{i},strcat(result,'emg.sto'))));
        ID_calib_AD(i) = importdata(strcat(path,strcat(trial{i},'ID.sto')));
        ID_calib_EMG(i) = importdata(strcat(path1,strcat(trial{i},'ID.sto')));
        ID_calib_AD(i).data(1,:) = [];
        ID_calib_EMG(i).data(1,:) = [];
        emgFilt_calib_AD(i) = importdata(strcat(path,strcat(trial{i},'emgFilt.sto')));
        emgFilt_calib_EMG(i) = importdata(strcat(path1,strcat(trial{i},'emgFilt.sto')));
    end

    JS_calib_AD = importdata(strcat(path_calib,'jointStiffness.sto'));
    JS_calib_EMG = importdata(strcat(path1_calib,'jointStiffness.sto'));

%% Segmenting torque of Calibration

    for i = 1:nRep
        ID_calib{i} = ID_calib_AD(i).data(:,3);
        ID_calib_1{i} = ID_calib_EMG(i).data(:,3);
        ID{i} = ID_AD(i).data(:,3);
        ID_1{i} = ID_EMG(i).data(:,3);
    end

    for i = 1:length(locs)
        for j = 1
            ID_calib_seg{i}{j} = ID_calib{i}(1:locs{i}(j),1);
            ID_calib_seg_1{i}{j} = ID_calib_1{i}(1:locs{i}(j),1);
        end
    end

    for i = 1:length(locs)
        for j = 2:length(locs{i})
            ID_calib_seg{i}{j} = ID_calib{i}(locs{i}(j-1)+1:locs{i}(j),1);
            ID_calib_seg_1{i}{j} = ID_calib_1{i}(locs{i}(j-1)+1:locs{i}(j),1);
        end
    end

    for i = 1:length(locs)
        ID_calib_seg{i}{1} = [];
        ID_calib_seg{i} = ID_calib_seg{i}(~cellfun('isempty',ID_calib_seg{i}));
        ID_calib_seg_1{i}{1} = [];
        ID_calib_seg_1{i} = ID_calib_seg_1{i}(~cellfun('isempty',ID_calib_seg_1{i}));
    end

    for i = 1:length(ID_calib_seg)
        for j = 1:length(ID_calib_seg{i})
            ID_calibx = length(ID_calib_seg{i}{j});
            ID_calib_pad{i}{j} = padarray(ID_calib_seg{i}{j},round(ID_calibx*0.3),'replicate','both');
            ID_calib_pad_1{i}{j} = padarray(ID_calib_seg_1{i}{j},round(ID_calibx*0.3),'replicate','both');
            ID_calib_1600{i}{j}  = resample(ID_calib_pad{i}{j},1600,length(ID_calib_pad{i}{j}));
            ID_calib_1600_1{i}{j}  = resample(ID_calib_pad_1{i}{j},1600,length(ID_calib_pad_1{i}{j}));
            ID_calib_1000{i}{j} = ID_calib_1600{i}{j}(301:1300);
            ID_calib_1000_1{i}{j} = ID_calib_1600_1{i}{j}(301:1300);
        end
    end

    for i = 1:nRep
        stiff = cell2mat(ID_calib_1000{i});
        ID_calib_1000_sub{i} = stiff;
        ID_calib_final{i} = mean(ID_calib_1000_sub{i},2);
        stiff_1 = cell2mat(ID_calib_1000_1{i});
        ID_calib_1000_sub_1{i} = stiff_1;
        ID_calib_final_1{i} = mean(ID_calib_1000_sub_1{i},2);
    end

%% Segmenting torque of CEINMS

    for i = 1:nRep
        ID{i} = ID_AD(i).data(:,3);
        ID_1{i} = ID_EMG(i).data(:,3);
    end

    for i = 1:length(locs)
        for j = 1
            ID_seg{i}{j} = ID{i}(1:locs{i}(j),1);
            ID_seg_1{i}{j} = ID_1{i}(1:locs{i}(j),1);
        end
    end

    for i = 1:length(locs)
        for j = 2:length(locs{i})
            ID_seg{i}{j} = ID{i}(locs{i}(j-1)+1:locs{i}(j),1);
            ID_seg_1{i}{j} = ID_1{i}(locs{i}(j-1)+1:locs{i}(j),1);
        end
    end

    for i = 1:length(locs)
        ID_seg{i}{1} = [];
        ID_seg{i} = ID_seg{i}(~cellfun('isempty',ID_seg{i}));
        ID_seg_1{i}{1} = [];
        ID_seg_1{i} = ID_seg_1{i}(~cellfun('isempty',ID_seg_1{i}));
    end

    for i = 1:length(ID_seg)
        for j = 1:length(ID_seg{i})
            idx = length(ID_seg{i}{j});
            ID_pad{i}{j} = padarray(ID_seg{i}{j},round(idx*0.3),'replicate','both');
            ID_pad_1{i}{j} = padarray(ID_seg_1{i}{j},round(idx*0.3),'replicate','both');
            ID_1600{i}{j}  = resample(ID_pad{i}{j},1600,length(ID_pad{i}{j}));
            ID_1600_1{i}{j}  = resample(ID_pad_1{i}{j},1600,length(ID_pad_1{i}{j}));
            ID_1000{i}{j} = ID_1600{i}{j}(301:1300);
            ID_1000_1{i}{j} = ID_1600_1{i}{j}(301:1300);
        end
    end

    for i = 1:nRep
        stiff = cell2mat(ID_1000{i});
        ID_1000_sub{i} = stiff;
        ID_final{i} = mean(ID_1000_sub{i},2);
        stiff_1 = cell2mat(ID_1000_1{i});
        ID_1000_sub_1{i} = stiff_1;
        ID_final_1{i} = mean(ID_1000_sub_1{i},2);
    end

%% Segmenting JS of CEINMS

    for i = 1:nRep
        JS{i} = JS_AD(i).data(:,3);
        JS_1{i} = JS_EMG(i).data(:,3);
    end

    for i = 1:length(locs)
        for j = 1
            JS_seg{i}{j} = JS{i}(1:locs{i}(j),1);
            JS_seg_1{i}{j} = JS_1{i}(1:locs{i}(j),1);
        end
    end

    for i = 1:length(locs)
        for j = 2:length(locs{i})
            JS_seg{i}{j} = JS{i}(locs{i}(j-1)+1:locs{i}(j),1);
            JS_seg_1{i}{j} = JS_1{i}(locs{i}(j-1)+1:locs{i}(j),1);
        end
    end

    for i = 1:length(locs)
        JS_seg{i}{1} = [];
        JS_seg{i} = JS_seg{i}(~cellfun('isempty',JS_seg{i}));
        JS_seg_1{i}{1} = [];
        JS_seg_1{i} = JS_seg_1{i}(~cellfun('isempty',JS_seg_1{i}));
    end

    for i = 1:length(JS_seg)
        for j = 1:length(JS_seg{i})
            idx = length(JS_seg{i}{j});
            JS_pad{i}{j} = padarray(JS_seg{i}{j},round(idx*0.3),'replicate','both');
            JS_pad_1{i}{j} = padarray(JS_seg_1{i}{j},round(idx*0.3),'replicate','both');
            JS_1600{i}{j}  = resample(JS_pad{i}{j},1600,length(JS_pad{i}{j}));
            JS_1600_1{i}{j}  = resample(JS_pad_1{i}{j},1600,length(JS_pad_1{i}{j}));
            JS_1000{i}{j} = JS_1600{i}{j}(301:1300);
            JS_1000_1{i}{j} = JS_1600_1{i}{j}(301:1300);
        end
    end

    for i = 1:nRep
        stiff = cell2mat(JS_1000{i});
        JS_1000_sub{i} = stiff;
        JS_final{i} = mean(JS_1000_sub{i},2);
        JS_SD{i} = std(JS_1000_sub{i},[],2);
        stiff_1 = cell2mat(JS_1000_1{i});
        JS_1000_sub_1{i} = stiff_1;
        JS_final_1{i} = mean(JS_1000_sub_1{i},2);
        JS_SD_1{i} = std(JS_1000_sub_1{i},[],2);
    end

%% R2 & RMSE of AD - Trial data

    for i = 1:nRep
        predicted = ID_AD(i).data(:,3);
        actual = ID_calib_AD(i).data(:,3);
        RMSE_ID_AD(i) = sqrt(sum(((predicted-actual).^2))/length(predicted));
        nRMSE_ID_AD(i) = RMSE_ID_AD(i)./rms(actual);
        R_ID_AD(i) = ((length(predicted).*sum(predicted.*actual))-(sum(predicted).*sum(actual)))/sqrt(((length(predicted).*sum(predicted.^2))-(sum(predicted).^2)).*((length(actual).*sum(actual.^2))-(sum(actual).^2)));
        R2_ID_AD(i) = (R_ID_AD(i)).^2;

        for j = 1:length(JS_1000{i})
            predicted = JS_1000{i}{j};
            actual = JS_calib_AD.data(:,3);
            RMSE_JS_AD{i}(j) = sqrt(sum(((predicted-actual).^2))/length(predicted));
            R_JS_AD{i}(j) = ((length(predicted).*sum(predicted.*actual))-(sum(predicted).*sum(actual)))/...
                sqrt(((length(predicted).*sum(predicted.^2))-(sum(predicted).^2)).*...
                ((length(actual).*sum(actual.^2))-(sum(actual).^2)));
            R2_JS_AD{i}(j) = (R_JS_AD{i}(j)).^2;
            nRMSE_JS_AD{i}(j) = RMSE_JS_AD{i}(j)/rms(actual);
        end
    end

%% RMSE & R2 of EMG - Trial data

    for i =1:nRep
        predicted = ID_EMG(i).data(:,3);
        actual = ID_calib_EMG(i).data(:,3);
        RMSE_ID_EMG(i) = sqrt(sum(((predicted-actual).^2))/length(predicted));
        R_ID_EMG(i) = ((length(predicted).*sum(predicted.*actual))-(sum(predicted).*sum(actual)))/sqrt(((length(predicted).*sum(predicted.^2))-(sum(predicted).^2)).*((length(actual).*sum(actual.^2))-(sum(actual).^2)));
        R2_ID_EMG(i) = (R_ID_EMG(i)).^2;
        nRMSE_ID_EMG(i) = RMSE_ID_EMG(i)/rms(actual);
        for j = 1:length(JS_1000_1{i})
            predicted = JS_1000_1{i}{j};
            actual = JS_calib_EMG.data(:,3);
            RMSE_JS_EMG{i}(j) = sqrt(sum(((predicted-actual).^2))/length(predicted));
            R_JS_EMG{i}(j) = ((length(predicted).*sum(predicted.*actual))-(sum(predicted).*sum(actual)))/sqrt(((length(predicted).*sum(predicted.^2))-(sum(predicted).^2)).*((length(actual).*sum(actual.^2))-(sum(actual).^2)));
            R2_JS_EMG{i}(j) = (R_JS_EMG{i}(j)).^2;
            nRMSE_JS_EMG{i}(j) = RMSE_JS_EMG{i}(j)/rms(actual);
        end
    end


%% Mean & SD of RMSE & R2 of Stiffness data

    for i = 1:nRep
        mean_nRMSE_JS_EMG(i) = mean(nRMSE_JS_EMG{i});
        mean_nRMSE_JS_AD(i) = mean(nRMSE_JS_AD{i});
        mean_R2_JS_EMG(i) = mean(R2_JS_EMG{i});
        mean_R2_JS_AD(i) = mean(R2_JS_AD{i});
        mean_RMSE_JS_EMG(i) = mean(RMSE_JS_EMG{i});
        mean_RMSE_JS_AD(i) = mean(RMSE_JS_AD{i});

        SD_nRMSE_JS_EMG(i) = std(nRMSE_JS_EMG{i});
        SD_nRMSE_JS_AD(i) = std(nRMSE_JS_AD{i});
        SD_R2_JS_EMG(i) = std(R2_JS_EMG{i});
        SD_R2_JS_AD(i) = std(R2_JS_AD{i});
        SD_RMSE_JS_EMG(i) = std(RMSE_JS_EMG{i});
        SD_RMSE_JS_AD(i) = std(RMSE_JS_AD{i});
    end

    R_sq_AD = cell2mat(R2_JS_AD);
    mean_R_sq = mean(R_sq_AD);
    SD_R_sq = std(R_sq_AD);
    RMSE_AD = cell2mat(RMSE_JS_AD);
    mean_RMSE = mean(RMSE_AD);
    SD_RMSE = std(RMSE_AD);
    nRMSE_AD = cell2mat(nRMSE_JS_AD);
    mean_nRMSE = mean(nRMSE_AD);
    SD_nRMSE = std(nRMSE_AD);
    R_sq_EMG = cell2mat(R2_JS_EMG);
    mean_R_sq = mean(R_sq_EMG);
    SD_R_sq = std(R_sq_EMG);
    RMSE_EMG = cell2mat(RMSE_JS_EMG);
    mean_RMSE = mean(RMSE_EMG);
    SD_RMSE = std(RMSE_EMG);
    nRMSE_EMG = cell2mat(nRMSE_JS_EMG);
    mean_nRMSE = mean(nRMSE_EMG);
    SD_nRMSE = std(nRMSE_EMG);

%% Plot for checking

   for i = 1:nRep
        figure
        plot(JS_final_1{i});
        hold on
        plot(JS_calib_EMG.data(:,3));
        legend('Stiffness-Estimated','Stiffness-Reference');
        xlabel('Index');
        ylabel('Stiffness (Nm/rad)');
    end

    for i = 1:nRep
        figure
        plot(ID_final_1{i});
        hold on
        plot(ID_calib_final_1{i});
        legend('Estimated Torque','Reference Torque');
            xlabel('time');
            ylabel('Torque');
    end

%% maximum & minimum stiffness

    for i = 1:nRep
        for j = 1:length(JS_1000{i})
            max_stiff_EMG = max(max(max(JS_1000_1{i}{j})));
            max_stiff_AD = max(max(max(JS_1000{i}{j})));
            min_stiff_EMG = min(min(min(JS_1000_1{i}{j})));
            min_stiff_AD = min(min(min(JS_1000{i}{j})));
            max_JS_calib_EMG = max(JS_calib_EMG.data(:,3));
            min_JS_calib_EMG = min(JS_calib_EMG.data(:,3));
        end
    end

%% Plot the matrices with error bars for their standard deviation

    % Calculate the standard deviation for each matrix
    sd1 = std(cell2mat(ID_final),[],2);
    sd3 = std(cell2mat(ID_calib_final),[],2);
    sd2 = std(cell2mat(ID_final_1),[],2);

    flex = linspace(0,100,1000);
    figure1 = figure;
        hold on;
        shadedErrorBar(flex, mean(cell2mat(ID_final),2), sd1, 'lineprops', '-b', 'patchSaturation', 0.2);
        shadedErrorBar(flex, mean(cell2mat(ID_final_1),2), sd2, 'lineprops', '-g', 'patchSaturation', 0.2);
        shadedErrorBar(flex, mean(cell2mat(ID_calib_final),2), sd3, 'lineprops', '-r', 'patchSaturation', 0.2);
        xlabel('flexion (%)');
        ylabel('Torque (Nm)');
        legend('Estimated-MU','Estimated-EMG','Reference');
        saveas(figure1,'Torque_final6.jpg')

    %% Calculate the standard deviation for each matrix

    % sd1 = std(cell2mat(JS_final_1),[],2);
    idx = length(stDev(:,10));
    sd_pad = padarray(stDev,round(idx*0.3),'replicate','both');   
    sd_1600  = resample(sd_pad,1600,length(sd_pad));
    sd_1000 = sd_1600(301:1300,:);

    sd6 = mean(sd_1000,2);
    sd4 = std(cell2mat(JS_SD),[],2);
    sd5 = std(cell2mat(JS_SD_1),[],2);

    %% Plot the matrices with shaded error bars for their standard deviation
    
    figure2 = figure;
        hold on;
        shadedErrorBar(flex, mean(cell2mat(JS_final),2), sd4, 'lineprops', '-b', 'patchSaturation', 0.2);
        shadedErrorBar(flex, mean(cell2mat(JS_final_1),2), sd5, 'lineprops', '-g', 'patchSaturation', 0.2);
        shadedErrorBar(flex, mean(JS_calib_AD.data(:,3),2), sd6, 'lineprops', '-r', 'patchSaturation', 0.2);
        xlabel('flexion (%)');
        ylabel('Stiffness (Nm/rad)');
        legend('Estimated-MU','Estimated-EMG','Reference');
        saveas(figure2,'Stiffness_final6.jpg')
    
    %% Save results
    
    stiff_calib = JS_calib_AD.data(:,3);
    torque_calib = ID_calib_final;
% 
    save('D:\SaiArunajatesan\Thesis\Data\subject01\Results_valid.mat','JS_final','JS_final_1','stiff_calib','ID_final','ID_final_1','torque_calib','sd1','sd2','sd3','sd4','sd5','sd6','R_sq_AD','RMSE_AD','nRMSE_AD','R_sq_EMG','RMSE_EMG','nRMSE_EMG','RMSE_ID_EMG','R2_ID_EMG','nRMSE_ID_EMG','RMSE_ID_AD','R2_ID_AD','nRMSE_ID_AD');