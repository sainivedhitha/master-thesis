%Create a figure with the plot of raw signals
%
%Input:
%   sig_path:       path of the file .sig
%   board_info:     structure storing the information 
%                   about the settings of the acquisition board
%                   Returned by load_xml.
%                   Refers to the help of load_xml for details
%   setup_info:     structure storing the information 
%                   about the setup used to acquire the signals
%                   Returned by load_xml.
%                   Refers to the help of load_xml for details
%   sig_info:       structure returned by load_xml.
%                   refers to the help of load_xml for details
%   s_start:        first second to read
%   epoch_length:   length of the signal epoch to load
%                   if not specified, the whole signal is loaded
%
%Output:
%
%   sig:            signal samples (in mV)
%                   cell array with one cell for each array used
%                   to acquire the signals
%                   Each cell contains the signals detected by an array
%                   (matrix NCHSxNSAMPS where
%                       NCHS: number of channels in the array
%                       NSAMPS: number of acquired samples) 
%
% Author: Marco Gazzoni
% Date: 20070319
%
%Last update:
%   20070320:   added the epoch_length input parameter
%               the signal epoch is loaded within this function
%
% Copyright (C) 27/10/2006 LISiN Politecnico di Torino (lisin@polito.it)
% All rights reserved

function sig= load_signal(sig_path, board_info, setup_info, sig_info, s_start, epoch_length)

%if the epoch length is not specified
if (nargin < 6)
    epoch_length= sig_info.length;
end

%load the signal
h = fopen([sig_path sig_info.filename]);
fseek(h, s_start*setup_info.nchs*board_info.fsamp*2, 'bof');
sig = fread(h, [setup_info.nchs board_info.fsamp*epoch_length], 'short');
fclose(h);
