function [shifts, minErrors, stDev] = alignment(f, fs, nr, nRep, ensemblesTorqDec, j_samples_dec, tau_samples_dec)

% Alignment algorithm aimed at reducing variability in sinusodial
% realizations
%
% Inputs 
% 
% f                 : task frequency [Hz]
% fs                : sampling frequency [Hz]
% nr                : Number of valid realizations per trial
% nRep              : Number of trials for specific condition
% ensembleTorqDec   : matrix containing torque ensemble (after decimation)
% j_samples_dec     : number of samples in averaging segments
% tau_samples_dec   : number specifying the maximum shift checked by the
%                   alignment algorithm
% 
% Outputs
% shifts            : matrix specifying the shifts (samples) calculated by the
%                   alignment algorithm (for each time point and realization)
% minErrors         : matrix specifying the minimum error calculated by the
%                   alignment algorithm (for each time point and realization)
% stDev             : matrix specifying the standard deviation of each
%                   realization
% 
% 
% Written by : Alejandro Moya Esteban
% Date       : May 2018


% Find parameters (amplitude and frequency) of the target sine signal
% Take only the actual period of the signal (not the extra added samples
% needed to shift window). 
y = ensemblesTorqDec(j_samples_dec+tau_samples_dec:ceil(fs/f)+j_samples_dec+tau_samples_dec-1,:);

target = mean(y,2);
ampl = max(target)-min(target);

i1 = find(target == max(target));
i2 = find(target == min(target));

period = 2*abs(i2-i1)/fs;
freq = 1/period;

% We use 1/f because we need a complete period of the signal
t = -j_samples_dec/fs:1/fs:1/f+j_samples_dec/fs;
sine = -ampl/2*sin(2*pi*freq*t+pi/2);
% plot(t,sine)

% Initialize vector/matrices needed for the alignment algorithm
shifts = zeros(ceil(fs/f),nRep*nr);                            % Matrix storing shifts for each realization and time point.         
minErrors = zeros(ceil(fs/f),nRep*nr);                              % Matrix needed to remove outliers based on MSE
stDev = zeros(ceil(fs/f),nRep*nr);                                  % Matrix needed to remove outliers based on STD
targetSine = zeros(2*j_samples_dec+1,ceil(fs/f));                   % Auxiliary variable to make the code more readable
realization = zeros(2*j_samples_dec+1,2*tau_samples_dec+1,nRep*nr); % Auxiliary variable to make the code more readable

% Alignment algorithm
'Alignment'
for i = 1:ceil(fs/f)-1 % Find shift for each time point
    i
    i_low = i;
    i_high = i+2*j_samples_dec;
    targetSine(:,i) = sine(i_low:i_high);
    for k = 1:nRep*nr     % Find shift for each realization in the ensemble
        error = zeros(2*tau_samples_dec+1,1);
        for tau = 1:2*tau_samples_dec+1     % Find the best shift over a range [-tau_samples_dec, tau_samples_dec]
            tau_low = tau+i-1;
            tau_high = tau+2*j_samples_dec+i-1;
            realization(:,tau,k) = ensemblesTorqDec(tau_low:tau_high,k);
            error(tau) = immse(targetSine(:,i), realization(:,tau,k));      % Compute MSE
        end
        shifts(i,k) = find(error == min(error));  % Find the shift for which the error is minimum                           
        if shifts(i,k) == 2*tau_samples_dec + 1
            shifts(i,k) = 2*tau_samples_dec;
        end
        minErrors(i,k) = min(error);              % Store the minimum error for specific time point and realization
        stDev(i,k) = std(ensemblesTorqDec(i+shifts(i,k):i+2*j_samples_dec+shifts(i,k),k));  % Store standard deviation for the aligned segment
    end
end

end