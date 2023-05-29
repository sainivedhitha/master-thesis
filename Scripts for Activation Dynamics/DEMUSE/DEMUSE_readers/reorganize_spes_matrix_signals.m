function [MONO SD]= reorganize_spes_matrix_signals(Signal)

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
IN_B= 2;
IN_C= 3;
IN_D= 4;
%   IN_A= 2;
%   IN_B= 1;
%   IN_C= 3;
%   IN_D= 4;

%calc the first channel (-1) the matrix connector is connected to
ch1_IN_A= (IN_A-1)*16;
ch1_IN_B= (IN_B-1)*16;
ch1_IN_C= (IN_C-1)*16;
ch1_IN_D= (IN_D-1)*16;

%order of the channel in the matrix connectors
o16c= ([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]);
o16cb= fliplr([8 7 6 5 4 3 2 1 16 15 14 13 12 11 10 9]);    %for flipped adapter
%o16cb= ([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]);

% %Electrode map top side view
%Electrode map top side view
Electrode_Map = ...
[ch1_IN_B+o16cb(9)   ch1_IN_A+o16c(1)    ch1_IN_A+o16c(2)    ch1_IN_A+o16c(3)    NaN;
ch1_IN_B+o16cb(11)   ch1_IN_B+o16cb(10)   ch1_IN_A+o16c(6)    ch1_IN_A+o16c(5)    ch1_IN_A+o16c(4); ...
ch1_IN_B+o16cb(13)   ch1_IN_B+o16cb(12)   ch1_IN_A+o16c(16)   ch1_IN_A+o16c(8)    ch1_IN_A+o16c(7); ...
ch1_IN_B+o16cb(8)    ch1_IN_B+o16cb(16)   ch1_IN_A+o16c(15)   ch1_IN_A+o16c(13)   ch1_IN_A+o16c(14); ...
ch1_IN_B+o16cb(4)    ch1_IN_B+o16cb(7)    ch1_IN_B+o16cb(14)   ch1_IN_A+o16c(12)   ch1_IN_A+o16c(11); ...
ch1_IN_B+o16cb(3)    ch1_IN_B+o16cb(6)    ch1_IN_B+o16cb(15)   ch1_IN_A+o16c(10)   ch1_IN_A+o16c(9); ...
ch1_IN_B+o16cb(1)    ch1_IN_B+o16cb(2)    ch1_IN_B+o16cb(5)    ch1_IN_C+o16c(2)    ch1_IN_C+o16c(1); ...
ch1_IN_D+o16c(11)   ch1_IN_D+o16c(10)   ch1_IN_D+o16c(9)    ch1_IN_C+o16c(3)    ch1_IN_C+o16c(4); ...
ch1_IN_D+o16c(12)   ch1_IN_D+o16c(13)   ch1_IN_C+o16c(8)    ch1_IN_C+o16c(6)    ch1_IN_C+o16c(5); ...
ch1_IN_D+o16c(15)   ch1_IN_D+o16c(14)   ch1_IN_C+o16c(14)   ch1_IN_C+o16c(16)   ch1_IN_C+o16c(7); ...
ch1_IN_D+o16c(7)    ch1_IN_D+o16c(8)    ch1_IN_D+o16c(16)   ch1_IN_C+o16c(13)   ch1_IN_C+o16c(15); ...
ch1_IN_D+o16c(4)    ch1_IN_D+o16c(5)    ch1_IN_D+o16c(6)    ch1_IN_C+o16c(11)   ch1_IN_C+o16c(12); ...
ch1_IN_D+o16c(3)    ch1_IN_D+o16c(2)    ch1_IN_D+o16c(1)    ch1_IN_C+o16c(9)    ch1_IN_C+o16c(10)]';

[n_rows,n_cols]=size(Electrode_Map);
for i_col=1:n_cols
    for i_row=1:n_rows
        if(isfinite(Electrode_Map(i_row,i_col)))
           MONO{i_row,i_col}=Signal(Electrode_Map(i_row,i_col),:);
        else
           MONO{i_row,i_col}= [];
        end
    end
end

for i_col=1:n_cols
    for i_row=1:n_rows-1
        if(isfinite(Electrode_Map(i_row+1,i_col)))
           SD{i_row,i_col}=Signal(Electrode_Map(i_row+1,i_col),:)-Signal(Electrode_Map(i_row,i_col),:);
        else
           SD{i_row,i_col}= [];
        end
    end
end


