function data = DEMUSEtool_reader_RAMPS_MIRIAM(filepath,filename,epoch_length)
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

% Sub = filename(4:5);
% Contract = filename(8);
% delimInd = findstr(filepath,'\');
% newFilePath = filepath(1:delimInd(end-3))
% 
% [rawsig, SIG, dd_mtx, lp_mtx, force]=Load_sig_hand_auto2(newFilePath,Sub,Contract,[1,75000]);
% clear rawsig dd_mtx lp_mtx;
% SIG=chooseChannelsCell(SIG,fsamp,1,0,50);

load([filepath filename]);

%f=fopen([filepath(1:end-length(filename)-1) filename]);
%[filepath(1:end-length(filename)-1) filename]
%EMG=fread(f,[64,inf],'short');

%fclose(f);

%eval(['EMG =',filename(1:end-4),';']); 

%EMG = EMG(find(good),:)';

%
for k1=1:13,
    data.SIG{k1,1} = EMG(k1,:)';
end
for k1=14:26,
    data.SIG{k1-13,2} = EMG(k1,:)';
end
for k1=27:39,
    data.SIG{k1-26,3} = EMG(k1,:)';
end
for k1=40:52,
    data.SIG{k1-39,4} = EMG(k1,:)';
end
for k1=53:64,
    data.SIG{k1-52,5} = EMG(k1,:)';
end

data.signal_length = size(EMG,2)
data.force =  [];
data.montage = 'SD';
data.IED = 10;
data.fsamp = 2048;
data.gain = 5000;
data.AUXchannels = [];