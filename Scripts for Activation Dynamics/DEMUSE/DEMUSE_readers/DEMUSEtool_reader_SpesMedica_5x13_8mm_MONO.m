function data = DEMUSEtool_reader_SpesMedica_5x13_8mm_MONO(filepath,filename,epoch_length)
% function data = DEMUSEtool_reader_SpesMedica_5x13_8mm_MONO(filepath,filename,epoch_length)
%
% DEMUSEtool reader: reads surface EMG acquired by a matrix of 5x13
% electrodes with inter-electrode distance of 8 mm in monopolar montage - for dynamic prtocols. 
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

if nargin < 3
    epoch_length = inf;
end

%obtain xml filename
xml_filename= [filename(1:end-6) '.xml'];
%load experimental session info
[session_info sbj_info board_info setup_info sig_info]= load_xml([filepath xml_filename]);

sig_num = str2num(filename(end-5:end-4));
data.SIG = load_signal(filepath, board_info, setup_info, sig_info{sig_num}, epoch_length);
data.SIG = data.SIG.*2500./2^11/sig_info{sig_num}.gain(1);
data.force = data.SIG(123,:);

data.SIG = reorganize_spes_matrix_signals2(data.SIG(1:64,:),'MONO');

data.signal_length = length(data.force);
data.montage = 'MONO';
data.IED = 8;
data.fsamp = board_info.fsamp;
%data.gain = 1;
data.AUXchannels = [];
data.AUXchannels_description = {};