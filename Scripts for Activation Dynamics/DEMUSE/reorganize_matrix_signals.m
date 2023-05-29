function [MONO, SD, DD]= reorganize_spes_matrix_signals(Signal, matrix)

% function MONO=reorganize_spes_matrix_signals(Signal)
%
% Input
%   Signal: Matrix that contains the signal 
%   matrix: type of matrix (SPES for spes medica matrix)
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

if(strcmpi(matrix, 'SPES'))
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
        ch1_IN_D+o16c(3)    ch1_IN_D+o16c(2)    ch1_IN_D+o16c(1)    ch1_IN_C+o16c(9)    ch1_IN_C+o16c(10)];
else
    %Electrode map top side view
    Electrode_Map = ...
        [NaN   3    2   1    NaN; ...
         4     5    6   7   8; ...
         13    12   11  10  9; ...
         14    15   16  17  18; ...
         23     22  21  20  19; ...
         24     25  26  27  28; ...
         33     32  31  30  29; ...
         34     35  36  37  38; ...
         43     42  41  40  39; ...
         44     45  46  47  48; ...
         53     52  51  50  49; ...
         54     55  56  57  58; ...
         nan    61  60  59  nan];
end

[n_rows,n_cols]=size(Electrode_Map);
MONO = [];

for i_col=1:n_cols
    for i_row=1:n_rows
        if(isfinite(Electrode_Map(i_row,i_col)))
           MONO{i_row,i_col}=Signal(Electrode_Map(i_row,i_col),:);
        else
           MONO{i_row,i_col}= [];
        end
    end
end

SD = [];
% for i_col=1:n_cols
%     for i_row=1:(n_rows-1)
%         if(isfinite(Electrode_Map(i_row,i_col)) && isfinite(Electrode_Map(i_row+1,i_col)))
%             SD{i_row,i_col} =  Signal(Electrode_Map(i_row,i_col),:) - Signal(Electrode_Map(i_row+1,i_col),:);          
%         else
%             SD{i_row,i_col}=[];
%         end
%     end
% end

DD= [];
% for i_col=1:n_cols
%     for i_row=1:(n_rows-2)
%         if(isfinite(Electrode_Map(i_row,i_col)))
%            DD{i_row,i_col}=SD{i_row,i_col}-SD{i_row+1,i_col};
%         else
%             DD{i_row,i_col}=[];
%         end
%     end
% end
