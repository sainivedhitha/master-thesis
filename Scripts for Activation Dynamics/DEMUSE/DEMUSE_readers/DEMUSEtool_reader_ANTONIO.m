function data = DEMUSEtool_reader_ANTONIO(filepath,filename,epoch_length)
% function data = DEMUSEtool_reader_iEMG_MultiChAA(filepath,filename,epoch_length)
%
% DEMUSEtool reader: reads surface EMG acquired by a multichannel intramuscular electrode
%                    developed at Aalborg University (SMI). 
%
% INPUTS:
%   - filepath: directory with the SIG file to be loaded
%   - finename: SIG file to be loaded
%   - epoch_length: (optional) length of the epoch of signal to be loaded (in s)
%
% OUTPUT:
%   data: structure with the following fields;
%       SIG - two dimensional cell array with surface EMG channel in each
%             cell - SIG{r,c} is the channel in row r and column c. Missing
%             electrodes are denoted by empty arrays, e.g. SIG{1,1} = [];
%       fsamp - sampling frequency of sEMG
%       signal_length - length of a surface EMG signals (in samples)
%       montage - montage of electrodes - 'MONO' for monopolar, 'SD' for
%                 single differential
%       IED - inter-electrode distance (in mm)
%       force - measured force signal if avalable, empty array otherwise
%       AUXchannels - auxilary channels (currently not used by DEMUSEtool
%       AUXchannels_description - cell array of texts describing the data
%                   in AUXchannels (one cell per channel)
%
% -----------------------------------------------------------------------
% Copyright: LISiN, Politecnico di Torino, Italy
%            SSL, FEECS, University of Maribor, Slovenia
% Author: Ales Holobar (ales.holobar@uni-mb.si)
% Last modified: 14. 10. 2008

fs = 2048;

D = load([filepath, filename]);
D = cell2mat(struct2cell(D));
Data = D.EMGData;
data.force = [];
data.gain=1000;
Data = Data(1:64, 5*fs:67*fs); % skipping the first 5 s and cutting the period after 67 s;
data.SIG = {};
count = 1;
Data = Data';
%% Notch filter 
for w = 50:50:500  % artifacts at these frequencies
    wo = w/(fs/2);  bw = wo/100;
    [b,a] = iirnotch(wo,bw);
    %fvtool(b,a);
    Data = filtfilt(b,a,Data);
end

Data = Data';

for c = 1: 8 %size(Ind,2)
    for r = 1 : 8 %size(Ind,1)
        if r < 9
            data.SIG{r,c} = Data(count,:);
            count = count + 1;
        else
            if count < 65
            data.SIG{r,c} = Data(count,:);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.signal_length = length(Data);
data.montage = 'MONO';
data.IED = 8;
data.fsamp = fs;
data.gain = 1000;
data.AUXchannels = [];
data.AUXchannels_description = {};
clear Data;