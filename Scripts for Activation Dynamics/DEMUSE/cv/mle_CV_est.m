function [cv,cvCorr] = mle_CV_est(Segna,dint,fsamp)

% Function for the estimation of CV (global or of single MU) using Maximum Likelihood
% multiple-channel methods. The computational cost is very low since the iterative Newton
% method is used for efficient minimum detection in a similar way as it is done for the 
% case of two channels in [2].
%
% Input parameters:
%
% Segna		matrix with the EMG signals in raws (no innervation zones)
% dint		interelectrode distance (in meters)
% fsamp		sampling frequency (in Hz)
%
% Output parameters:
%
% cv			the estimated conduction velocity 
%
% References:
%
% [1] D. Farina, W. Muhammad, E. Fortunato, O. Meste, R. Merletti, H. Rix, "Estimation of 
% single motor unit conduction velocity from the surface EMG signal detected 
% with linear electrode arrays", Medical & Biological Eng. & Comput
% (development of the technique)
%
% [2] Mc Gill K,  Dorfman L. High resolution alignment of sampled waveforms. IEEE Trans Biomed Eng 1984;31:462-70.
% (Maximum likelihood estimator for only two channels) 
%
% Method developed by D. Farina, E. Fortunato, R. Merletti, O. Meste, W. Muhammad, H. Rix
% Software written by D. Farina and W. Muhammad
% Reference e-mail: farina@athena.polito.it

% Flip of the signal matrix in order to have always a correct direction of propagation
primo = 1; secondo = 2;
if size(Segna,1)>3
primo = 2; secondo = 3;
end

corr=xcorr(Segna(primo,:),Segna(secondo,:));

[a b]=max(corr);
if b > ceil(length(corr)/2)
   Segna=flipud(Segna);
end;

% Estimate the starting point for the iterative process. Inportant not to fix
% an absolute starting point (see [1] for details)
start = localac(Segna(primo,:),Segna(secondo,:),dint,fsamp)
% Application of the iterative procedure to estimate CV
[cv,teta] = mle3(Segna,start,dint,fsamp);

% added by Ales - estimation of channel correlation coefficients 
interpFact=10;
iSenga(1,:) = [zeros(1,100*interpFact) interp(Segna(1,:),interpFact) zeros(1,100*interpFact)];  
for k1=2:size(Segna,1);
  iSenga(k1,:) = [zeros(1,(100-round(k1*teta))*interpFact) interp(Segna(k1,:),interpFact) zeros(1,(100+round(k1*teta))*interpFact)];    
end
tmpCC1 = corrcoef(iSenga');
for k1=2:size(Segna,1);
  iSenga(k1,:) = [zeros(1,(100+round(k1*teta))*interpFact) interp(Segna(k1,:),interpFact) zeros(1,(100-round(k1*teta))*interpFact)];    
end
tmpCC2 = corrcoef(iSenga');
ind = find(tmpCC1<1);
cvCorr = median(max([tmpCC1(ind)'; tmpCC2(ind)']));
%cv=cv*sign(tmpCC2(1,2)-tmpCC1(1,2));
% END of added by Ales