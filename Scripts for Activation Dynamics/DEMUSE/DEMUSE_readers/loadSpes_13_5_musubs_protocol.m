%function to view the signals acquired during a session
%using the EMG-USB amplifier and the Acquisition software
%
%   signal_path: path of the directory that contains the signals acquired
%               during the experimental session
%   session_code:   code of the experimental session
%   isig:   ordinal number of the signal acquired during the experimental
%           session the user want to check
%   sstart: first second to read
%   epoch_length: length of the epoch of signal to load (in s)
%
%ex:
%   MONO= check_signal_quality('..\..\signals\dynamic-elbow-flex-ext-01\measures\dynamic-elbow-flex-ext-01\', '08MaGaLi0705100929', 3, 30);
%function [SIG, ElbowAngle, fsamp, gain]= loadPermutatedSignal2(signal_path, session_code, isig, epoch_length)
function [MONO, SD, fsamp, gain]= loadSpes_13_5_musubs_protocol(signal_path, session_code, isig, s_start, epoch_length)

if (nargin < 5)
    epoch_length= 1;
end

%obtain xml filename
xml_filename= [session_code '.xml'];

%load experimental session info
[session_info sbj_info board_info setup_info sig_info]= load_xml([signal_path xml_filename]);

nsigs= length(sig_info);

sig= load_signal([signal_path], board_info, setup_info, sig_info{isig}, s_start, epoch_length);
% conversion to mV
sig = sig.*2500./2^11;
sig = sig/sig_info{isig}.gain(1);



% [Be,Ae] = butter (4, [20/board_info.fsamp*2 350/board_info.fsamp*2]); % EMG band pass filter
% for ch=1:size(sig,1)
%     sig(ch,:)= filtfilt(Be,Ae,sig(ch,:));
% end

%Reorganize the acquired channel in a cell array {row,col}
%each cell contains the channel samples
[MONO SD]= reorganize_spes_matrix_signals(sig(1:64,:));

% %Reorganize the acquired channel in a cell array {row,col}
% %each cell contains the channel samples
% SIG{1,5} = [];
% c=5; for r=2:13
%     SIG{r,c} = sig(r-1,:);
% end
% c=4; for r=1:13
%     SIG{14-r,c} = sig(12+r,:);
% end
% c=3; for r=1:13
%     SIG{r,c} = sig(25+r,:);
% end
% c=2; for r=1:13
%     SIG{14-r,c} = sig(38+r,:);
% end
% c=1; for r=1:13
%     SIG{r,c} = sig(51+r,:);
% end
% SIG= SIG';
%[n_rows,n_cols]=size(SD);

% for i_col=1:n_cols
%     for i_row=1:n_rows-1
%         if(~isempty(SIG{i_row+1,i_col}))
%            SD{i_row,i_col}=SIG{i_row+1,i_col}-SIG{i_row,i_col};
%         else
%            SD{i_row,i_col}= [];
%         end
%     end
% end

fsamp = board_info.fsamp;
gain = sig_info{isig}.gain(1);