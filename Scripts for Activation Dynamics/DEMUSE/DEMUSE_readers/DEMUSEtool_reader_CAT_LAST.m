function data = DEMUSEtool_reader_CAT_LAST(filepath,filename,epoch_length)
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

load([filepath filename(1:end-4)]);

f = fopen([filepath filename]);
EMG=fread(f,[NCHANNEL,inf],'short');
fclose(f);

[BF,AF]=butter(4,[1000 4900]*2/10240);
iEMG = filtfilt(BF,AF,EMG(1:8,:)')';

AUX = EMG(73:end,1:10240/2048:end);

[BF,AF]=butter(4,[10 500]*2/10240);
EMG(9:72,:)=filtfilt(BF,AF,EMG(9:72,:)')';
EMG=EMG(9:72,1:10240/2048:end);
sEMG=EMG;



Finterf =60;

Fvibration = 138;

Fsamp=2048;
sigLength=length(EMG);
Frad=round(Fsamp/200);
fInt=round(0:Finterf*sigLength/Fsamp:(sigLength/2));
fInt2=round(0:Fvibration*sigLength/Fsamp:(sigLength/2));

c=1;
for k1=13:-1:1,
    EMG(k1,:)=real(subtractFreq2(EMG(k1,:),fInt,Frad,1));
    EMG(k1,:)=real(subtractFreq2(EMG(k1,:),fInt2,Frad,0));
    pause
    data.SIG{c,1} = EMG(k1,:);
    index = find(EMG(k1,1:end-1)>2000 & EMG(k1,2:end)<2000);
    if length(index)>2, [c,1],data.SIG{c,1}=[];end;
    c=c+1;
end
c=1;
for k1=14:1:26,
    EMG(k1,:)=real(subtractFreq2(EMG(k1,:),fInt,Frad,0));
    EMG(k1,:)=real(subtractFreq2(EMG(k1,:),fInt2,Frad,0));
    data.SIG{c,2} = EMG(k1,:)';
    index = find(EMG(k1,1:end-1)>2000 & EMG(k1,2:end)<2000);
    if length(index)>2, [c,1],data.SIG{c,2}=[];end;
    c=c+1;
end
c=1;
for k1=39:-1:27,
    EMG(k1,:)=real(subtractFreq2(EMG(k1,:),fInt,Frad,0));
    EMG(k1,:)=real(subtractFreq2(EMG(k1,:),fInt2,Frad,0));
    data.SIG{c,3} = EMG(k1,:)';
    index = find(EMG(k1,1:end-1)>2000 & EMG(k1,2:end)<2000);
    if length(index)>2, [c,1],data.SIG{c,3}=[];end;    
    c=c+1;
end
c=1;
for k1=40:52,
    EMG(k1,:)=real(subtractFreq2(EMG(k1,:),fInt,Frad,0));
    EMG(k1,:)=real(subtractFreq2(EMG(k1,:),fInt2,Frad,0));
    data.SIG{c,4} = EMG(k1,:)';
    index = find(EMG(k1,1:end-1)>2000 & EMG(k1,2:end)<2000);
    if length(index)>2, [c,1],data.SIG{c,4}=[];end;
    c=c+1;
end
c=1;
for k1=53:64,
    EMG(k1,:)=real(subtractFreq2(EMG(k1,:),fInt,Frad,0));
    EMG(k1,:)=real(subtractFreq2(EMG(k1,:),fInt2,Frad,0));
    data.SIG{c,5} = EMG(k1,:)';
    index = find(EMG(k1,1:end-1)>2000 & EMG(k1,2:end)<2000);
    if length(index)>2, [c,1],data.SIG{c,5}=[];end;
    c=c+1;
end

data.signal_length = size(EMG,2);
data.force =  [];
data.montage = 'SD';
data.IED = 10;
data.fsamp = 2048;
data.gain = 5000;
data.AUXchannels = [];

save([filepath filename(1:end-4),'_raw'],'iEMG','sEMG','AUX','NCHANNEL','COMMENT');

clear iEMG sEMG AUX NCHANNEL COMMENT EMG;