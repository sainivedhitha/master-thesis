function data = DEMUSEtool_reader_SpesMedica_5x13_8mm_SD_permutator_64chs_oldXML(filepath,filename,epoch_length)
% function data = DEMUSEtool_reader_SpesMedica_5x13_8mm_SD_permutator_64chs_oldXML(filepath,filename,epoch_length)
%
% DEMUSEtool reader: reads surface EMG acquired by a matrix of 5x13
% electrodes with inter-electrode distance of 8 mm - configuration from Aalborg
% with old XML file.
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

[SIG, force] = dynAutoLoad64ch(filepath,filename,1,2048);
%SIG=chooseChannelsCell(SIG,fsamp,1,0,50);

data.SIG = SIG;
data.signal_length = length(force);
data.force = force;
data.montage = 'SD';
data.IED = 8;
data.fsamp = 2048;
%data.gain = 1;
data.AUXchannels = [];
data.AUXchannels_description = {};