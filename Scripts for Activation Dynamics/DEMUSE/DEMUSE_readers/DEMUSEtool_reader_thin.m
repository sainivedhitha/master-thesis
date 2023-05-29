function data = DEMUSEtool_reader_thin(filepath,filename,epoch_length)
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

f=fopen([filepath filename]);
EMG=fread(f,[64,inf],'short');
fclose(f);

clc

%[massimo,index] = max(s,[],2);

%find(massimo>2000)

%figure,plot(s(massimo>2000,:)');

%if ~isempty(massimo>2000),disp('SATURATION!!!!!!!');end

EMG=double(EMG)*5/2^12;

%EMG=EMG(1:16,:);
%
c=1;
for k1=1:16,    
    data.SIG{c,1} = EMG(k1,:);
    if max(EMG(k1,:))>2000, [c,1],data.SIG{c,1}=[];end;
    c=c+1;
end

data.signal_length = size(EMG,2);
data.force =  EMG(25,:);
data.montage = 'SD';
data.IED = 10;
data.fsamp = 10000;
data.gain = 5000;
data.AUXchannels = [];

clear EMG;