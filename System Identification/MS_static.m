% function [ensemblesPosDec, ensemblesTorqDec, ensemblesPertDec, shifts, minErrors, stDev, useRealization, segmentLength,...
%     HI, K, stiffness, j_samples_dec, tau_samples_dec] = MS_static(order, subject, bstrapMax)

% Function which based on alignment and MS algorithm estimates joint
% stiffness in static situation (fixed ankle position but TV torque)
%
% Inputs 
% 
% order             : number specifying the order in which the subject 
%                   performed the static task. 1st (1), 2nd (2) or 3rd (3)
% subject           : string specifying the subject of study, e.g. 'Subject5'
% btstrpMax         : maximum number of iteration of the MS algorithm
%
% Outputs
% 
% ensemblesPosDec   : matrix containing position ensemble (after decimation)
% ensembleTorqDec   : matrix containing torque ensemble (after decimation)
% ensemblePertDec   : matrix containing perturbation ensemble (after
%                   decimation)
% shifts            : matrix specifying the shifts calculated by the
%                   alignment algorithm (for each time point and realization)
% minErrors         : matrix specifying the minimum error calculated by the
%                   alignment algorithm (for each time point and realization)
% stDev             : matrix specifying the standard deviation of each
%                   realization
% useRealization    : matrix with values either 0 (realization is an
%                   outlier) or 1 (realization is OK)
% segmentLength     : vector containing segment lengths used by the MS
%                   algorithm (indicated in seconds)
% HI                : cell structure containing double sided IRF's
% K                 : cell structure containing stiffness estimates
% stiffness         : matrix containing stiffness estimates (averaged)
% j_samples_dec     : number of samples in averaging segments
% tau_samples_dec   : number specifying the maximum shift checked by the
%                   alignment algorithm
% 
% Written by : Alejandro Moya Esteban
% Date       : May 2018

% Definition of parameters 
clear all
close all
clc

nRep = 8;                                           % Number of trials for specific condition
fs = 2048;                                          % Sampling frequency [Hz]
dt = 1/fs;                                          % Sampling time [s]
f = 0.5;                                            % Sine frequency [Hz]
df = 20;                                            % Decimation factor
nlags = floor(0.12*fs/df);                          % Number of lags considered in the IRF estimation
t1 = 43;                                            % Initial time point (avoid trans. effects)
t2 = 77;                                           % Final time point (avoid final exp. effects)
samplesTrial = (t2-t1)*fs+1;                        % Total number of samples in a trial
nr = (t2-t1)*f-2;                                   % Number of valid realizations per trial
jmax = 10;                                          % Choose the max segment length
extraSamples = 3*nlags+jmax*5/2;                    % 3*nlags is a requirement of the MS algorithm

%%
% Ludvig files should be added to the path in order to run the MS
% algorithm. Modify the following line according to current location of the
% files. 
addpath(pwd,'D:\SaiArunajatesan\Thesis\Multisegment algorithm - Stiffness via system identification\DL_Tools')

% Specify the path in which raw experimental data is located. 
directory = 'D:\SaiArunajatesan\Thesis\Data\subject06\';

% Since the order of the conditions was randomized across subjects, the
% variable order is necessary to specify the position (1st, 2nd or 3rd) in
% which the specific condition was performed. 
% if order == 1
    trial = 'Trial3C1R';
    trial1 = 'PRBS_6-';    
% elseif order == 2
%     trial = '\D02_OfflineTVIdent\Trial10C1R';
% else
%     trial = '\D02_OfflineTVIdent\Trial13C1R';
% end

% Generate the ensembles with exact number of samples
j_samples = ceil(0.25*fs)+3*nlags*df;

% If signals are decimated, 1 second can be used because there is not a
% huge effect in the computation time. However, note that very few
% realizations (if any) will need a 1-second shift. 
tau_samples = ceil(1*fs);

% Initialization of ensembles
ensemblesTorq = zeros(ceil(fs/f)+(j_samples+tau_samples)*2,nRep*nr);
ensemblesPos = zeros(ceil(fs/f)+(j_samples+tau_samples)*2,nRep*nr);
ensemblesPert = zeros(ceil(fs/f)+(j_samples+tau_samples)*2,nRep*nr);
% Initialization of torque vector
torque = zeros(samplesTrial,nRep);

%%
% Segmentation of recordings. Each iteration of the 'for' loop corresponds to
% one 2-minute trial in the experimental session. 
for i = 1:nRep
    load(strcat(directory, strcat(trial, num2str(i))));
    Perturb = load(strcat('D:\SaiArunajatesan\Thesis\Data\subject06\', strcat(trial1, num2str(i))));

    % Load measured position, measured torque and perturbation signal
    posComplete = OutputData.AchData(:,2);
    torqueComplete = OutputData.AchData(:,4);
    perturbationComplete = Perturb.PRBS_trick;
    
    % Load measured torque (removing the first ant last 3 seconds)
    torque(:,i) = OutputData.AchData(t1*fs:t2*fs,4);
    
    % Delay between the position perturbation and the actual position
    % (needed to correctly estimate stiffness)
    maxLag = 1200;                                 % In samples
    Phiyu = xcorr(posComplete, perturbationComplete, maxLag);
    delay = find(Phiyu == max(Phiyu))-maxLag;

    % Determine torque minimum peaks
    minPeakHeight = 2.5;
    [pks,locs] = findpeaks(-torque(:,i), 'MinPeakHeight', minPeakHeight,'MinPeakDistance', 0.7*ceil(fs/f));
    segments = diff(locs);
    
    % Create ensembles and resample all realizations so that all present
    % the same number of samples. 
    x = linspace(0,100,ceil(fs/f)+(j_samples+tau_samples)*2);
    aux = (j_samples + tau_samples);               % Define auxiliary variable in order to make the code more readable
    for j =1:nr
        ensemblesPos(:,(i-1)*nr+j) = spline(linspace(0,100,segments(j)+aux*2+1),posComplete(locs(j)+t1*fs-aux:locs(j+1)+t1*fs+aux), x);
        ensemblesTorq(:,(i-1)*nr+j) = spline(linspace(0,100,segments(j)+aux*2+1),torqueComplete(locs(j)+t1*fs-aux:locs(j+1)+t1*fs+aux), x);
        ensemblesPert(:,(i-1)*nr+j) = spline(linspace(0,100,segments(j)+aux*2+1),perturbationComplete(locs(j)+t1*fs-aux-delay:locs(j+1)+t1*fs+aux-delay), x);
    end
end

%%

[r1,c1] = size(ensemblesPert);
figure,
plot(ensemblesPert(:,1:c1));

figure,
plot(ensemblesTorq(:,1:c1));

figure,
plot(ensemblesPos(:,1:c1));

%%
% Decimation of realizations
fs = fs/df;
j_samples_dec = ceil(0.25*fs)+3*nlags;
tau_samples_dec = ceil(1*fs)-1;

% Initialize ensembles (for decimated signals)
ensemblesTorqDec = zeros(ceil(length(ensemblesTorq)/df),nRep*nr);
ensemblesPosDec = zeros(ceil(length(ensemblesPos)/df),nRep*nr);
ensemblesPertDec = zeros(ceil(length(ensemblesPert)/df),nRep*nr);

% Perform the actual decimation on each ensemble
for i = 1:nRep*nr
    ensemblesTorqDec(:,i) = decimate(ensemblesTorq(:,i),df);
    ensemblesPosDec(:,i) = decimate(ensemblesPos(:,i),df);
    ensemblesPertDec(:,i) = decimate(ensemblesPert(:,i),df);
end

%%

[r1,c1] = size(ensemblesPertDec);

figure,
plot(ensemblesPertDec(:,1:c1));

figure,
plot(ensemblesTorqDec(:,1:c1));

figure,
plot(ensemblesPosDec(:,1:c1));

%%
% Alignment algorithm
[shifts, minErrors, stDev] = alignment(f, fs, nr, nRep, ensemblesTorqDec, j_samples_dec, tau_samples_dec);

% Outlier removal
useRealization = outlierRemoval(f, fs, nr, nRep, minErrors, stDev);

%% Multisegment algorithm
'Multisegment algorithm: bootstrap'
bstrapMax = 35;
K = cell(bstrapMax,jmax);

for bstrp = 1:bstrapMax             % Repeat the MS algorithm over bstrapMax iterations
    strcat('BOOTSTRAP ITERATION NUMBER: ', num2str(bstrp))
    for i = 1:ceil(fs/f)-1            % One stiffness estimate per time point.
        i   % To keep track of the algorithm in the command window. 
        % Initialize position, torque and perturbation matrices based on
        % the number of outliers at each time point. 
        nr_outliers = sum(useRealization,2);
        posd = zeros(2*j_samples_dec+1,nr_outliers(i));
        tqd = zeros(2*j_samples_dec+1,nr_outliers(i));
        ptd = zeros(2*j_samples_dec+1,nr_outliers(i));
        index = 0;
        for k = 1:nr*nRep
            if useRealization(i,k) == 1     % Avoid using outliers.
                index = index+1;
                lag = shifts(i,k);
                
                % Select aligned realizations
                posd(:,index) = ensemblesPosDec(i+lag:i+2*j_samples_dec+lag,k);
                tqd(:,index) = ensemblesTorqDec(i+lag:i+2*j_samples_dec+lag,k);
                ptd(:,index) = ensemblesPertDec(i+lag:i+2*j_samples_dec+lag,k);
            end
        end            

        for j = 1:jmax   % Repeat the estimation varying the window length, which is determined by 'j'
            rs = randsample(nr_outliers(i),ceil(0.95*nr_outliers(i)));  % Choose randomly 95 % of the realizations.
            l1 = extraSamples-floor(j*5/2);
            l2 = extraSamples+floor(j*5/2);
            
            % Ludvig's multisegment algorithm. Note that the mean of the
            % torque and position ensembles is removed to highlight
            % effects induced by the PRBS perturabation signal. 
            [Hi, Yi, vaft] = multistiffcl(posd(:,rs)-mean(posd(:,rs),2),tqd(:,rs)-mean(tqd,2),ptd(:,rs),l1+1,l2+1,1/fs,1,nlags);
            HI{j}(:,i) = Hi;
        end  
    end
    for j = 1:jmax
        K{j}(:,bstrp) = sum(HI{j}).*(1/fs);     % Integrate the double sided IRF's to estimate stiffness. 
    end
end

%% Matrix containing stiffness estimates for each segment length
stiffness = zeros(ceil(fs/f)-1, jmax);
for j =1:jmax
    stiffness(:,j) = -mean(K{j},2);        % Average stiffness over 35 iterations of the MS algorithm.
end

% Vector with the length of each averaging segment.
segmentLength = zeros(1,jmax);
for j = 1:jmax
    l1 = floor(j*5/2);
    segmentLength(j) = (2*l1+1)/fs;
end

% end

%%

[r1,c1] = size(posd);

figure,
plot(posd(:,1:c1));

figure,
plot(tqd(:,1:c1));

figure,
plot(ptd(:,1:c1));
%% save data - static
% 
path = 'D:\SaiArunajatesan\Thesis\Data\subject06\';
name = 'Static';
save(strcat(path,strcat(name,'.mat')),'bstrapMax', 'ensemblesPosDec', 'ensemblesTorqDec', 'ensemblesPertDec', 'shifts', 'minErrors', 'stDev', 'useRealization', 'segmentLength', 'HI', 'K', 'stiffness', 'j_samples_dec', 'tau_samples_dec');
