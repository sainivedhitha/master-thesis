%% Script to run simulation studies
close all
clear all
clc

% Load position PRBS perturbation signal
% load('PRBS.mat')

trial1 = 'PRBS_6-';   
fs = 2048;
dt = 1/fs;
ids = [round(41*fs),round(80*fs)];

% Load torque and prepare it for the simulation

for i = 8
    PRBS_file = load(strcat('D:\SaiArunajatesan\Thesis\Data\subject06\', strcat(trial1, num2str(i))));
    PRBS = PRBS_file.PRBS_trick;
    t = 0:1/fs:length(PRBS)/fs-1/fs;
    Tmax = max(t);
    torque = 5*sin(2*pi*0.5*t)'; % Theoretical sine
end
    % Create time series inputs for the Simulink model. 
    torque_TS = timeseries(torque,t);
    PRBS_TS = timeseries(PRBS,t);

    % Virtual parameters of the admittance controller of the Achilles
    I = 0.1;
    D = 2.5;
    K = 4500;

    A = [I D K];
    B = 1;

    % Create PRBS_trick (filter & inverse of the admittance controller)
    H = tf(B,A);
    H_z = c2d(H, dt);
    H_inv = H^-1;
    H_inv_z = c2d(H_inv, dt, 'tustin');

    fc = 40;
    [D,C] = maxflat(0,2,2*fc*dt);
    H_filter_z = tf(D,C,dt);

    H_final_z = H_filter_z*H_inv_z;

    B_filt = H_final_z.Numerator;
    B_filt = B_filt{1};
    A_filt = H_final_z.Denominator;
    A_filt = A_filt{1};

    PRBS_trick = filter(B_filt, A_filt, PRBS);
    figure,
    subplot(211)
    plot(t,PRBS)
    xlim([40 80])
    subplot(212)
    plot(t,PRBS_trick)
    xlim([40 80])
    
    sim('simulationModel')
    set_param('simulationModel/Admittance Controller','Denominator','[A(1) A(2) A(3)]');

    % Retrieve data from Simulink
    measTorque = Output_AC.Data + 0.2*randn(size(Output_AC.Data));
    pos_pert = pos_pert.Data;
    pos_noPert = pos_noPert.Data;
    sim_stiffness = sim_stiffness.Data;

    figure(3)
    subplot(311)
    plot(t,measTorque)
    xlim([40 80])
    title('Measured torque')

    subplot(312)
    plot(t, pos_noPert)
    xlim([40 80])
    title('Non-perturbed position (admittance position)')

    subplot(313)
    plot(t, pos_pert)
    xlim([40 80])
    title('Position + perturbation')

    figure(4)
    subplot(211)
    plot(t,measTorque)
    xlim([40 80])
    title('Measured torque')
    subplot(212)
    plot(t,sim_stiffness)
    xlim([40 80])
    title('Simulated Stiffness')


    %%
    % Now include a multisegment algorithm to check the estimated stiffness
    addpath(pwd,'D:\SaiArunajatesan\Thesis\Multisegment algorithm - Stiffness via system identification\DL_Tools')

    % Important parameters
    f = 0.5;                            % Signal frequency [Hz.]
    fs = 2048;                          % Achilles' sampling frequency [Hz.]
    df = 20;                            % Decimation factor
    nlags = floor(0.12*fs/df);          % Number of lags considered in the IRF estimation

    % Define the segment of the signal consider for analysis
    t1 = 10;                            % Lower limit [s.]
    t2 = 118;                           % Upper limit [s.]
    nr = (t2-t1)*f;                     % Number of periods/realiz in the interval

    % Find delay (in samples) between position and PRBS perturbation
    maxLag = 1200;                     % In samples
    [Phiyu, tau] = xcorr(pos_pert, PRBS, maxLag);
    delay = find(Phiyu == max(Phiyu))-maxLag;

    % Generate the ensembles with exact number of samples
    j_samples = ceil(0.25*fs)+3*nlags*df;

    % If signals are decimated, 1 second can be used because there is not a
    % huge effect in the computation time. However, note that very few
    % realizations (if any) will need a 1-second shift. 
    tau_samples = ceil(1*fs);

    periods = t1:1/f:t2;
    % The 1.25 comes from 1 s. due to tau1 and 0.25 due to J1
    ensemblesTorq = zeros(ceil(fs/f)+(j_samples+tau_samples)*2,nr);
    ensemblesPos = zeros(ceil(fs/f)+(j_samples+tau_samples)*2,nr);
    ensemblesPert = zeros(ceil(fs/f)+(j_samples+tau_samples)*2,nr);
    for j = 1:nr
        low = periods(j);
        high = periods(j+1);
        ensemblesTorq(:,j) = measTorque(fs*low-(j_samples+tau_samples):fs*high+j_samples+tau_samples-1); 
        ensemblesPos(:,j) = pos_pert(fs*low-(j_samples+tau_samples):fs*high+j_samples+tau_samples-1);
        ensemblesPert(:,j) = PRBS(fs*low-(j_samples+tau_samples)-delay:fs*high+j_samples+tau_samples-delay-1);
    end

    fs = fs/df;
    j_samples_dec = ceil(0.25*fs)+3*nlags;
    tau_samples_dec = ceil(1*fs)-1;

    ensemblesTorqDec = zeros(ceil(length(ensemblesTorq)/df),nr);
    ensemblesPosDec = zeros(ceil(length(ensemblesPos)/df),nr);
    ensemblesPertDec = zeros(ceil(length(ensemblesPert)/df),nr);
    for i = 1:nr
        ensemblesTorqDec(:,i) = decimate(ensemblesTorq(:,i),df);
        ensemblesPosDec(:,i) = decimate(ensemblesPos(:,i),df);
        ensemblesPertDec(:,i) = decimate(ensemblesPert(:,i),df);
    end

    % Find parameters (amplitude and frequency) of the target sine signal
    % Take only the actual period of the signal (not the extra added samples
    % needed to shift windows. 
    y = ensemblesTorqDec(j_samples_dec+tau_samples_dec:ceil(fs/f)+j_samples_dec+tau_samples_dec-1,:);
    % plot(mean(y,2))

    target = mean(y,2);
    ampl = max(target)-min(target);

    i1 = find(target == max(target));
    i2 = find(target == min(target));

    period = 2*abs(i2-i1)/fs;
    freq = 1/period;

    % We use 1/f because we need a complete period of the signal (0.8 Hz),
    % (that's the time, t, that we will estimate) although in this time vector 
    % a sine with a different frequency is going to be inserted.

    t = -j_samples_dec/fs:1/fs:1/f+j_samples_dec/fs;
    sine = -ampl/2*sin(2*pi*freq*t+pi/2);
    % plot(t,sine)

    shifts = zeros(ceil(fs/f),nr);
    targetSine = zeros(2*j_samples_dec+1,ceil(fs/f));
    realization = zeros(2*j_samples_dec+1,2*tau_samples_dec+1,nr);

    'Alignment'
    for i = 1:ceil(fs/f)-1
        i
        i_low = i;
        i_high = i+2*j_samples_dec;
        targetSine(:,i) = sine(i_low:i_high);
        for k = 1:nr
            error = zeros(2*tau_samples_dec+1,1);
            for tau = 1:2*tau_samples_dec+1
                tau_low = tau+i-1;
                tau_high = tau+2*j_samples_dec+i-1;
                realization(:,tau,k) = ensemblesTorqDec(tau_low:tau_high,k);
                error(tau) = immse(targetSine(:,i), realization(:,tau,k));
            end
            shifts(i,k) = find(error == min(error));
        end
    end

    %%
    % shifts = shifts - tau_samples;
    aux = shifts - tau_samples_dec;
    timeShift = aux/fs;

    jmax = 10;                              % Choose the max segment length
    extraSamples = 3*nlags+jmax*5/2;             % The 36 is 3*nlags.

    % Just to have a vector with the length of each averaging segment.
    segmentLength = zeros(1,jmax);
    for j = 1:jmax
        l1 = floor(j*5/2);
        segmentLength(j) = (2*l1+1)/(fs/20);
    end

    'Multisegment algorithm'
    for i = 1:ceil(fs/f)-1
        i
        posd = zeros(2*j_samples_dec+1,nr);
        tqd = zeros(2*j_samples_dec+1,nr);
        ptd = zeros(2*j_samples_dec+1,nr);
        for k = 1:nr
            lag = shifts(i,k);
            lag = tau_samples_dec; % This involves 0 shift.
            posd(:,k) = ensemblesPosDec(i+lag:i+2*j_samples_dec+lag,k);
            tqd(:,k) = ensemblesTorqDec(i+lag:i+2*j_samples_dec+lag,k);
            ptd(:,k) = ensemblesPertDec(i+lag:i+2*j_samples_dec+lag,k);
        end
        for j = 1:jmax
            l1 = extraSamples-floor(j*5/2);
            l2 = extraSamples+floor(j*5/2);
    %         [Hi, Yi, vaft] = multistiffcl(posd-mean(posd,2),tqd-mean(tqd,2),ptd,l1+1,l2+1,1/fs,1,nlags);
            [Hi, Yi, vaft] = multistiffcl(posd-mean(posd,2),tqd-mean(tqd,2),ptd,l1+1,l2+1,1/fs,1,nlags);
            HI{j}(:,i) = Hi;
            VAF(j) = vaft;
        end  
    end
    K = cell(1,jmax);
    % With only one round (k) of the algorithm, it is not necessary to have
    % the cell structure. However, it is kept like this in case more rounds
    % of the algorithm want to be run (only if random realizations are
    % chosen instead of choosing them all). 
    k = 1;
    for j = 1:jmax
        K{j}(:,k) = sum(HI{j})*(1/fs);
    end

    stiffness = zeros(ceil(fs/f)-1, jmax);
    for j =1:jmax
        stiffness(:,j) = -mean(K{j},2);
    end

dir=strcat(strcat('D:\SaiArunajatesan\Thesis\Multisegment algorithm - Stiffness via system identification\Simulation files\Simulation Results\trial_6-8','.mat'));
save(dir,'torque_TS','PRBS_TS','PRBS_trick','measTorque','sim_stiffness','pos_noPert','pos_pert','ensemblesTorqDec','ensemblesPosDec','ensemblesPertDec','period','sine','realization','segmentLength','posd','tqd','ptd','K','stiffness');
