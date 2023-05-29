function useRealization = outlierRemoval(f, fs, nr, nRep, minErrors, stDev)
% Outlier removal based on min. errors and standard deviation

% Inputs 
% 
% f                 : task frequency [Hz]
% fs                : sampling frequency [Hz]
% nr                : Number of valid realizations per trial
% nRep              : Number of trials for specific condition
% minErrors         : matrix specifying the minimum error calculated by the
%                   alignment algorithm (for each time point and realization)
% stDev             : matrix specifying the standard deviation of each
%                   realization
% Output
% useRealization    : matrix with values either 0 (realization is an
%                   outlier) or 1 (realization is OK)

Y = interp1(linspace(1/numel(minErrors),1,length(minErrors)), sort(minErrors), 0.9, 'nearest');
Y = prctile(minErrors,90,2);
stdLow = prctile(stDev,5,2);
stdHigh = prctile(stDev,95,2);

useRealization = zeros(ceil(fs/f),nRep*nr);          % 0 = outlier, 1 = OK
for i = 1:ceil(fs/f)
    for k = 1:nRep*nr
        if minErrors(i,k) < Y(i) && stDev(i,k) > stdLow(i) && stDev(i,k) < stdHigh(i)
            useRealization(i,k) = 1;
        end
    end
end

end