function SIG = reorganize_spes_matrix_signals2(Signal,El_Configuration)

% function MONO=reorganize_spes_matrix_signals(Signal)
%
% Input
%   Signal: Matrix that contains the signal 
% Output
%   MONO = cell array that contains the reorganized monopolar signals
%
% Function that reorganize the monopolar signals acquired with 
%the Spes medica matrix
%
%Date: 20-4-2007
%Author: Troiano Amedeo
%Modified by Marco Gazzoni
%Last update: 20070509

%EMG-USB connector to which the 4 matrix connectors (A, B, C, D)
%are connected to
IN_A= 1;
IN_B= 4;
IN_C= 2;
IN_D= 3;
%calc the first channel (-1) the matrix connector is connected to
ch1_IN_A= (IN_A-1)*16;
ch1_IN_B= (IN_B-1)*16;
ch1_IN_C= (IN_C-1)*16;
ch1_IN_D= (IN_D-1)*16;

%order of the channel in the matrix connectors
o16c= [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];

%Electrode map top side view
Electrode_Map = ...
[ch1_IN_B+o16c(9)   ch1_IN_A+o16c(1)    ch1_IN_A+o16c(2)    ch1_IN_A+o16c(3)    NaN;
ch1_IN_B+o16c(11)   ch1_IN_B+o16c(10)   ch1_IN_A+o16c(6)    ch1_IN_A+o16c(5)    ch1_IN_A+o16c(4); ...
ch1_IN_B+o16c(13)   ch1_IN_B+o16c(12)   ch1_IN_A+o16c(16)   ch1_IN_A+o16c(8)    ch1_IN_A+o16c(7); ...
ch1_IN_B+o16c(8)    ch1_IN_B+o16c(16)   ch1_IN_A+o16c(15)   ch1_IN_A+o16c(13)   ch1_IN_A+o16c(14); ...
ch1_IN_B+o16c(4)    ch1_IN_B+o16c(7)    ch1_IN_B+o16c(14)   ch1_IN_A+o16c(12)   ch1_IN_A+o16c(11); ...
ch1_IN_B+o16c(3)    ch1_IN_B+o16c(6)    ch1_IN_B+o16c(15)   ch1_IN_A+o16c(10)   ch1_IN_A+o16c(9); ...
ch1_IN_B+o16c(1)    ch1_IN_B+o16c(2)    ch1_IN_B+o16c(5)    ch1_IN_C+o16c(2)    ch1_IN_C+o16c(1); ...
ch1_IN_D+o16c(11)   ch1_IN_D+o16c(10)   ch1_IN_D+o16c(9)    ch1_IN_C+o16c(3)    ch1_IN_C+o16c(4); ...
ch1_IN_D+o16c(12)   ch1_IN_D+o16c(13)   ch1_IN_C+o16c(8)    ch1_IN_C+o16c(6)    ch1_IN_C+o16c(5); ...
ch1_IN_D+o16c(15)   ch1_IN_D+o16c(14)   ch1_IN_C+o16c(14)   ch1_IN_C+o16c(16)   ch1_IN_C+o16c(7); ...
ch1_IN_D+o16c(7)    ch1_IN_D+o16c(8)    ch1_IN_D+o16c(16)   ch1_IN_C+o16c(13)   ch1_IN_C+o16c(15); ...
ch1_IN_D+o16c(4)    ch1_IN_D+o16c(5)    ch1_IN_D+o16c(6)    ch1_IN_C+o16c(11)   ch1_IN_C+o16c(12); ...
ch1_IN_D+o16c(3)    ch1_IN_D+o16c(2)    ch1_IN_D+o16c(1)    ch1_IN_C+o16c(9)    ch1_IN_C+o16c(10)]

[n_rows,n_cols]=size(Electrode_Map);
if strcmp(upper(El_Configuration),'MONO')    
    for i_col=1:n_cols
        for i_row=1:n_rows
            if(isfinite(Electrode_Map(i_row,i_col)))
                SIG{i_row,i_col}=Signal(Electrode_Map(i_row,i_col),:);
            else
                SIG{i_row,i_col}= [];
            end
        end
    end
elseif strcmp(upper(El_Configuration),'SD')
    for i_col=1:n_cols
        for i_row=1:(n_rows-1)
            if(isfinite(Electrode_Map(i_row,i_col)))                
                SIG{i_row,i_col}=Signal(Electrode_Map(i_row,i_col),:) - Signal(Electrode_Map(i_row+1,i_col),:);
            else
                SIG{i_row,i_col}=[];
            end
        end
    end    
elseif strcmp(upper(El_Configuration),'DD')
    for i_col=1:n_cols
        for i_row=1:(n_rows-2)
            if(isfinite(Electrode_Map(i_row,i_col)))                
                SIG{i_row,i_col}=Signal(Electrode_Map(i_row,i_col),:) - 2*Signal(Electrode_Map(i_row+1,i_col),:) + Signal(Electrode_Map(i_row+2,i_col),:);
            else
                SIG{i_row,i_col}=[];
            end
        end
    end
elseif strcmp(upper(El_Configuration),'LP')
    for i_col=1 : n_cols-2
        for i_row=1 : n_rows-2
            if(isfinite(Electrode_Map(i_row,i_col)))                
                SIG{i_row,i_col}=Signal(Electrode_Map(i_row,i_col+1),:) - 4*Signal(Electrode_Map(i_row+1,i_col+1),:) ...
                    + Signal(Electrode_Map(i_row+2,i_col+1),:) + Signal(Electrode_Map(i_row+1,i_col),:) ...
                    + Signal(Electrode_Map(i_row+1,i_col+2),:);
            else
                SIG{i_row,i_col}=[];
            end
        end
    end
end