function Y=chooseChannelsCell(y,Fsamp,cut,draw,Finterf)
%
% reformats the signals from data structure, supported by the LISiN load_sig method, to the signal matrix,
% supported by the CKC method. With parameter cut > 0 it also filters out the interfrences at frequency Finterf 
% and all its heigher harmonics.
% INPUTS:
%   y - input signal (data structure supported by load_sig method);
%   Fsamp - sampling frequency [Hz]
%   cut - set cut>0  to filter out the interferences at frequency Finterf
%   draw - set draw>0 to display the frequency content of original and filtered data
%   Finterf - basic artefact frequency [Hz] (e.g. for line interference Finterf=50 Hz)
% OUTPUT
%   Y - (filtered) signal matrix, with each signal in separate row

Frad=round(Fsamp/200);

h=waitbar(0,'Extracting and filtering the channels...');
count = 0;
sigLength = 0;

% Calculate number and the maximal length of channels

for k=1:size(y,1)        
    for m=1:size(y,2)                
        if (~isempty(y{k,m}))
            sigLength=max(sigLength,length(y{k,m}));
        end
    end
end

for k=1:size(y,1)        
    for m=1:size(y,2)                
        if (~isempty(y{k,m}))
            Y{k,m} = zeros(1,sigLength);
        end
    end
end

if nargin<5
    Finterf = 50;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter the channels
fInt=round(0:Finterf*sigLength/Fsamp:(sigLength/2));

if cut>0        
    for k=1:size(y,1)
        for m=1:size(y,2)                              
            if (~isempty(y{k,m}))                   
                ynew=real(subtractFreq2(y{k,m},fInt,Frad,draw));                   
                count = count +1;
                Y{k,m} = ynew;                       
            end
            waitbar(((k-1)*size(y,2)+m)/(size(y,1)*size(y,2)),h);
        end
    end
else
    for k=1:size(y,1)
        for m=1:size(y,2) 
            if (~isempty(y{k,m}))
                count = count + 1;
                Y(count,:) = y{k,m};        
            end            
            waitbar(((k-1)*size(y,2)+m)/(size(y,1)*size(y,2)),h);
        end
    end
end
close(h);