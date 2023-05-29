function data = DEMUSEtool_reader_MASS_UCCI_GASmedial(filepath,filename,epoch_length)
% function data = DEMUSEtool_reader_SpesMedica_5x13_8mm_SD(filepath,filename,epoch_length)
%
% DEMUSEtool reader: reads surface EMG acquired by a matrix of 5x13
% electrodes with inter-electrode distance of 8 mm. 
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
fsamp = 2048;

h = fopen([filepath filename]);
SIG = fread(h, [258 fsamp*epoch_length], 'short');
fclose(h);

data.force = SIG(257,:);


% xml_filename = [filename(1:end-6) '.xml']; %obtain xml filename
% [session_info sbj_info board_info setup_info sig_info]= load_xml([filepath xml_filename]); %load experimental session info
% sig_num = str2num(filename(end-5:end-4));
% try
%     gain = sig_info{sig_num}.gain(1);
% catch
%     gain = 1000;
% end
 
gain=1000;
DINAMICA = 10;			%	-5V +5V
PRECISION = 12; 		%	number of bits of the A/D converter
SIG = SIG.*DINAMICA/((2^PRECISION)*gain).*1e6; % The data are now converted in microvolt

%data.SIG{1,1} = [];


%%% arange channels into matrix data.SIG
%Ind = [ 1:12; 25:-1:14; 27:38; 51:-1:40; 53:64]';
%Ind = [ 1,1:11; 24:-1:13; 26:37; 50:-1:39; 52:63]';
data.SIG = {};
count = 193;
for r= 1: 13 %size(Ind,2)
    for c= 1 : 5 %size(Ind,1)
            data.SIG{r,c} = SIG(count,:);
            if count>256,
                data.SIG{r,c} = [];
            end
            count = count + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear SIG;
%data.SIG{1,1} = [];

data.signal_length = length(data.force);
data.montage = 'MONO';
data.IED = 8;
data.fsamp = fsamp;
%data.gain = 1;
data.AUXchannels = [];
%data.AUXchannels_description = sig_info{sig_num}.comments;