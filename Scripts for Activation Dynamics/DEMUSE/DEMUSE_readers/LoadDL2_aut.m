function [sig] = LoadDL2_aut(GAIN,reffile,FileLim);

if nargin<3
    FileLim=inf;
end

DINAMICA = 10;			%	-5V +5V
PRECISION = 12; 		%	number of bits of the A/D converter

% File loading
% [nomefile,path] = uigetfile('*.sig', 'Open datalogger SIG file');
% if nomefile == 0,  return;  end;
% reffile = [path nomefile];

reffile
h = fopen(reffile, 'r');
sig = fread(h, [16,FileLim], 'short');

% Arrange dynamic
sig = sig.*DINAMICA/((2^PRECISION)*GAIN).*1e6; % The data are now converted in microvolt