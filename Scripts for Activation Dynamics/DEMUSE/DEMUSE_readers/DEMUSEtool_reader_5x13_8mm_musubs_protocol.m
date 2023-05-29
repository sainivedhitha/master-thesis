function data = DEMUSEtool_reader_5x13_8mm_musubs_protoco(filepath,filename,epoch_length)
% function data = DEMUSEtool_reader_5x13_8mm_musubs_protocol(filepath,filename,epoch_length)
%
% DEMUSEtool reader: reads surface EMG acquired by a matrix of 5x13
% electrodes with inter-electrode distance of 8 mm - configuration for 
% dynamic surface EMG with elbow joint angle measured ba goniometer. 
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
% Adapted by Marco Gazzoni
% Last modified: 26.04.2011

%obtain xml filename
xml_filename= [filename(1:end-6) '.xml'];
%load experimental session info
[session_info sbj_info board_info setup_info sig_info]= load_xml([filepath xml_filename]);
[MONO, SD, fsamp, EMG_gain]= loadSpes_13_5_musubs_protocol(filepath, xml_filename(1:end-4), str2num((filename(end-5:end-4))), 0, 60)
%[SIG, ElbowAngle, fsamp, EMG_gain] = loadPermutatedSignal2(SIGFilePath,SIGFileName(1:end-6),str2num((SIGFileName(end-5:end-4))),[]);

data.SIG = SD;
data.signal_length = length(data.SIG{2,2});
data.force = [];
data.montage = 'SD';
data.IED = 8;
data.fsamp = fsamp;
%data.gain = EMG_gain;
data.AUXchannels = [];
data.AUXchannels_description = {};