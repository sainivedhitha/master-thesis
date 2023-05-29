function [sig, sdl_mtx, dd_mtx, lp_mtx, Force]=Load_sig_auto(path,subjectName,contraction,epoch);
% loads the signals from the dynamic protocol:
% IN:
%  path - path to the directory with sig files (e.g. C:\signals\)
%  subjectName - name acronyms of the subject whose SIG files are to be
%               loaded (e.g. 'AH');
%  contraction - character array denoting the contraction type (e.g. '5', '10', '30', 'I1', 'I2', 'I3', 'D1', 'D2', 'D3')
%  epoch - character array describing the epoch number (e.g. '1' '2' '3' '4' '5' '6' '7' '8' for contraction types '5' and '10'
%           '' (empty string) for contraction types '30', 'I1', 'I2', 'I3', 'D1', 'D2' and 'D3')
%
% OUT:
%   sig - raw signals (one signal in each raw)
%   sdl_mtx - cell structure with recorded longitudinal single differential signals  sdl_mtx{r,c}
%            corresponds to channel in row r and column c
%   dd_mtx - cell structure with recorded longitudinal double differential signals  dd_mtx{r,c}
%            corresponds to channel in row r and column c
%   lp_mtx - cell structure with recorded "Laplacian" signals lp_mtx{r,c}
%            corresponds to channel in row r and column c
%   Force - recorded force

%close all
warning off
%clc

%FSamp = repeat_answer('input(''FSamp ? '')');
FSamp = 2500;

%GAIN = repeat_answer('input(''GAIN ? '')');
GAIN=5000;

if ~isempty(epoch)
    epoch=['_' epoch];
end

% Built the string with all array's signal
[sig] = LoadDL2_aut(GAIN,[path 'PC1\' subjectName contraction epoch '_PC1.SIG']);		sig1 = sig;
% ...to remove the noise
for ch = 1 : size(sig,1) - 1
    sig1(ch,:) = remove_int(sig1(ch,:), FSamp, 15, 400);
end

[sig] = LoadDL2_aut(GAIN,[path 'PC2\' subjectName contraction epoch '_PC2.SIG']);		sig2 = sig;
% ...to remove the noise
for ch = 1 : size(sig,1) - 1
    sig2(ch,:) = remove_int(sig2(ch,:), FSamp, 15, 400);
end

[sig] = LoadDL2_aut(GAIN,[path 'PC3\' subjectName contraction epoch '_PC3.SIG']);		sig3 = sig;
% ...to remove the noise
for ch = 1 : size(sig,1) - 1
    sig3(ch,:) = remove_int(sig3(ch,:), FSamp, 15, 400);
end

[sig] = LoadDL2_aut(GAIN,[path 'PC4\' subjectName contraction epoch '_PC4.SIG']);		sig4 = sig;
% ...to remove the noise
for ch = 1 : size(sig,1) - 1
    sig4(ch,:) = remove_int(sig4(ch,:), FSamp, 15, 400);
end
Force=sig4(16,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           PRIMA PARTE                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resystem dimension signal
[sig, vpp, epocs] = adjst2(sig1, sig2, sig3, sig4,FSamp);

clear sig1 sig2 sig3 sig4;

% Buit the spatial filter (only LSD)
[sig_sdl] = filt_spat(sig);

% Transform mat to matcell
sdl_mtx = Resystem2(sig_sdl);

clear sig_sdl;
% Buit the spatial filter (the other)
dd_mtx = filt_spat2(sdl_mtx);
%[sdt_mtx, ddt_mtx, dd_mtx, ndd_mtx, ib2_mtx] = filt_spat2(sig, FSamp);

lp_mtx = filt_spat3(sdl_mtx,sig);

%clear sig

% Plot LSD
%plot_sig5(FSamp, vpp);

% Plot (the other)
%plot_sig6(sdt_mtx, ddt_mtx, dd_mtx, ndd_mtx, ib2_mtx, FSamp, vpp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       SECONDA PARTE                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LSD
% MNF and ARV
%global sdl_mtx;
%[arv_cell, fmean_cell] = param2_cell(sdl_mtx, FSamp, epocs);
% Build vector
%[arv_vett, fmean_vett] = build_vett_parameter(arv_cell, fmean_cell, epocs);
% Calc regression line about MNF and ARV parameter
%regression(arv_vett, fmean_vett, epocs);

% LDD
% MNF and ARV
%[arv_cell, fmean_cell] = param2_cell(dd_mtx, FSamp, epocs);
% Build vector
%[arv_vett, fmean_vett] = build_vett_parameter(arv_cell, fmean_cell, epocs);
% Calc regression line about MNF and ARV parameter
%regression(arv_vett, fmean_vett, epocs);

% NDD
% MNF and ARV
%[arv_cell, fmean_cell] = param2_cell_full(ndd_mtx, FSamp, epocs);
% Build vector
%[arv_vett, fmean_vett] = build_vett_parameter(arv_cell, fmean_cell, epocs);
% Calc regression line about MNF and ARV parameter
%regression(arv_vett, fmean_vett, epocs);

% IB2
% MNF and ARV
%[arv_cell, fmean_cell] = param2_cell(ib2_mtx, FSamp, epocs);
% Build vector
%[arv_vett, fmean_vett] = build_vett_parameter(arv_cell, fmean_cell, epocs);
% Calc regression line about MNF and ARV parameter
%regression(arv_vett, fmean_vett, epocs);