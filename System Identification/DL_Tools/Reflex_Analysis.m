%Code for finding reflexes in PRBS data

% clear old entries
clear all

% Set threshold for detecting perturbation to half the size of the
% perturbation
threshold = 0.03;

% Create vector of trials numbers
trials = [4:13];

for j = 1:length(trials)
    
    %Load trials
    trialno = trials(j);
    eval(['load Trial' num2str(trialno) '.mat'])
    
    %Transform data if necessary
    M = ([0 0 -1;-1 0 0;0 1 0]*[Mx,My,Mz]')';
    
    %Filter Tq at 125 Hz if necessary
    [b,a] = butter(4,1/100);
    Mf = filtfilt(b,a,M(:,3));
    
    posM = Position;
    tqM = Mf;
    
    % Remove mean and rectify EMG
    brdM = abs(detrend(EMG1,'constant'));
    bicM = abs(detrend(EMG2,'constant'));
    triLoM = abs(detrend(EMG3,'constant'));
    triLatM = abs(detrend(EMG8,'constant'));
   
    %Find positive going perturbations
    [ri,rj] = find((posM(1:end-1,:) < threshold) & (posM(2:end,:) > threshold));
    for i = 1:length(ri)
        if ri(i) > 250 & ri(i) < length(Position)-750
            % Gather position, torque and EMG for each perturbation
            posU(:,i) = posM(ri(i)-250:ri(i)+750,rj(i));
            tqU(:,i) = tqM(ri(i)-250:ri(i)+750,rj(i));
            bicU(:,i) = bicM(ri(i)-250:ri(i)+750,rj(i));
            brdU(:,i) = brdM(ri(i)-250:ri(i)+750,rj(i));
            triLoU(:,i) = triLoM(ri(i)-250:ri(i)+750,rj(i));
            triLatU(:,i) = triLatM(ri(i)-250:ri(i)+750,rj(i));
        end
    end
    
    %Find negative going perturbations
    [ri,rj] = find((posM(1:end-1,:) > threshold) & (posM(2:end,:) < threshold));
    for i = 1:length(ri)
       if ri(i) > 250 & ri(i) < length(Position)-750
            posD(:,i) = posM(ri(i)-250:ri(i)+750,rj(i));
            tqD(:,i) = tqM(ri(i)-250:ri(i)+750,rj(i));
            bicD(:,i) = bicM(ri(i)-250:ri(i)+750,rj(i));
            brdD(:,i) = brdM(ri(i)-250:ri(i)+750,rj(i));
            triLoD(:,i) = triLoM(ri(i)-250:ri(i)+750,rj(i));
            triLatD(:,i) = triLatM(ri(i)-250:ri(i)+750,rj(i));
        end
    end
    
    %Store mean of perturbation responses
    posUM(:,j) = mean(posU,2);
    tqUM(:,j) = mean(tqU,2);
    brdUM(:,j) = mean(brdU,2);
    triLoUM(:,j) = mean(triLoU,2);
    triLatUM(:,j) = mean(triLatU,2);
    bicUM(:,j) = mean(bicU,2);
    
    posDM(:,j) = mean(posD,2);
    tqDM(:,j) = mean(tqD,2);
    brdDM(:,j) = mean(brdD,2);
    triLoDM(:,j) = mean(triLoD,2);
    triLatDM(:,j) = mean(triLatD,2);
    bicDM(:,j) = mean(bicD,2);
end

%save results
%save PRBSresults posDM posUM tqDM tqUM bicDM bicUM brdDM brdUM triLoDM triLoUM triLatDM triLatUM 