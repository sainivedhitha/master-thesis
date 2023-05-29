% This script is to open Monopolar 64ch EMG matrix data.

function [Sig_OK, AUXsig] = auto_open_64_monopolar(n_ch,File_name,File_path)  %Change the number of measured ch in n_ch) (currently in total 124ch)

% Copyright (C) 19/07/2006 LISiN Politecnico di Torino (lisin@polito.it)
% All rights reserved
% Author: Alberto Botter (alberto.botter@delen.polito.it)
%
% Software for the organization of EMG signals detected in SD configuration
% with a 13*5 electrode array (Spes Medica) (13 rows aligned with muscle fibers)
% connected to the EMG128 in this way (skin side is UP!):
%       
%       |      * * * *    |
%       |    * * * * *    |
%       |    * * * * *    |                 A -> Array In 1
%   A ==|    * * * * *    |== B             B -> Array In 2 
%       |    * * * * *    |                 C -> Array In 3
%       |    * * * * *    |                 D -> Array In 4
%       |    * * * * *    |
%       |    * * * * *    |
%       |    * * * * *    |
%       |    * * * * *    |
%   C ==|    * * * * *    |== D
%       |    * * * * *    |
%       |    * * * * *    |

%
% Input: n_ch -> number of channels (EMG and AUX) saved in *.sig file
% Output: Sig_OK -> cell array. The cell "n" contains all the SD
%                   channels (from 1 to 12) of the column "n"
% 
%       Column n:   1 2 3 4 5
%
%       SD ch 1    NaN* * * *
%       SD ch 2     * * * * *
%       SD ch 3     * * * * * 
%       SD ch 4     * * * * *
%       SD ch 5     * * * * *
%       SD ch 6     * * * * *
%       SD ch 7     * * * * *
%       SD ch 8     * * * * *
%       SD ch 9     * * * * *
%       SD ch 10    * * * * *
%       SD ch 11    * * * * *
%       SD ch 12    * * * * *

% Select EMG signal

%% Open the file
h = fopen([File_path File_name],'r');
Signal = fread(h,[n_ch,inf],'short');
fclose all;

% sampling frequency
fsamp = 2048;

% Connectors configuration (see *.ppt files)
Ocm = [1 2 3 4];  %'New'/'current' Configuration
% Ocm = [2 1 4 3];   %'Old' Configuration

Electrode_Map =[(Ocm(3)-1)*16+10 (Ocm(3)-1)*16+09 (Ocm(4)-1)*16+01 (Ocm(4)-1)*16+02 (Ocm(4)-1)*16+03;...
    (Ocm(3)-1)*16+12 (Ocm(3)-1)*16+11 (Ocm(4)-1)*16+06 (Ocm(4)-1)*16+05 (Ocm(4)-1)*16+04;...
    (Ocm(3)-1)*16+15 (Ocm(3)-1)*16+13 (Ocm(4)-1)*16+16 (Ocm(4)-1)*16+08 (Ocm(4)-1)*16+07;...
    (Ocm(3)-1)*16+07 (Ocm(3)-1)*16+16 (Ocm(3)-1)*16+14 (Ocm(4)-1)*16+14 (Ocm(4)-1)*16+15;...
    (Ocm(3)-1)*16+05 (Ocm(3)-1)*16+06 (Ocm(3)-1)*16+08 (Ocm(4)-1)*16+13 (Ocm(4)-1)*16+12;...
    (Ocm(3)-1)*16+04 (Ocm(3)-1)*16+03 (Ocm(4)-1)*16+09 (Ocm(4)-1)*16+10 (Ocm(4)-1)*16+11;...
    (Ocm(3)-1)*16+01 (Ocm(3)-1)*16+02 (Ocm(2)-1)*16+05 (Ocm(2)-1)*16+02 (Ocm(2)-1)*16+01;...
    (Ocm(1)-1)*16+09 (Ocm(1)-1)*16+10 (Ocm(2)-1)*16+15 (Ocm(2)-1)*16+06 (Ocm(2)-1)*16+03;...
    (Ocm(1)-1)*16+11 (Ocm(1)-1)*16+12 (Ocm(2)-1)*16+14 (Ocm(2)-1)*16+07 (Ocm(2)-1)*16+04;...
    (Ocm(1)-1)*16+14 (Ocm(1)-1)*16+13 (Ocm(1)-1)*16+15 (Ocm(2)-1)*16+16 (Ocm(2)-1)*16+08;...
    (Ocm(1)-1)*16+07 (Ocm(1)-1)*16+08 (Ocm(1)-1)*16+16 (Ocm(2)-1)*16+12 (Ocm(2)-1)*16+13;...
    (Ocm(1)-1)*16+04 (Ocm(1)-1)*16+05 (Ocm(1)-1)*16+06 (Ocm(2)-1)*16+10 (Ocm(2)-1)*16+11;...
    (Ocm(1)-1)*16+nan (Ocm(1)-1)*16+03 (Ocm(1)-1)*16+02 (Ocm(1)-1)*16+01 (Ocm(2)-1)*16+09];

Electrode_Map = flipud(Electrode_Map)';

sig_dur = floor(size(Signal,2)/fsamp);
clear Sig_OK

% reconstruction of the channels
for columna = 1:size(Electrode_Map,1) % for each column
    for rowa = 1:size(Electrode_Map,2)  % for each row
        Sig_OK{rowa,columna} = zeros(1,length(Signal));
        if sum(Electrode_Map(columna,rowa)) < inf
            Sig_OK{rowa,columna} = Signal(Electrode_Map(columna,rowa),:);           
        else
            Sig_OK{rowa,columna} = zeros(1,length(Signal));
        end
     end
end

% AUX-channels 
if n_ch==124
    AUXsig(1,:) = Signal(121,:); %Force Sync
    AUXsig(2,:) = Signal(122,:); %TMS OUT
    AUXsig(3,:) = Signal(123,:); %Elect Stim OUT
    AUXsig(4,:) = Signal(124,:); %Force
elseif n_ch==128
    AUXsig(1,:) = Signal(122,:); %Force Sync
    AUXsig(2,:) = Signal(121,:); %Force KIN
    AUXsig(3,:) = Signal(124,:); %Elect Stim OUT
    AUXsig(4,:) = Signal(128,:); %Force
end
% To get mV: Divide y-value with amplification. The range is +- 2047(8)
    
    
    
    
    
    
    

