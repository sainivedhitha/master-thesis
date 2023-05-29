function [sig_out, force] = autoLoad64ch(path,filename,GAIN,FSamp,Hlim)
% loads the signals from the dynamic protocol:
% IN:
%  path - path to the directory with sig files (e.g. C:\signals\)
%  filename - name of the sig file to be loaded (without the path to this file, e.g. ramp_dec20060608180716.sig)
%  GAIN - signal gain (see experimental forms)
%  Fsamp - sampling frequency (2048 in this case)
%  Hlim - (optional) the length of loaded signals [is seconds], e.g, Hlim =
%           10 loads just the first 10 seconds of the recorded signals
%
% OUT:
%   sig_out - cell structure with recorded signals  sig_out{r,c}
%            corresponds to channel in row r and column c
%   force - recorded force

DINAMICA = 5;			%	-5V +5V
PRECISION = 12; 		%	number of bits of the A/D converter

if nargin<5
    Hlim=inf;
end

% File loading
%[filename,path] = uigetfile('*.sig', 'Open datalogger SIG file');
%if filename == 0,  return;  end;

reffile = [path filename]
h = fopen(reffile, 'r');
sig = fread(h, [64, Hlim], 'short');
fclose(h);

% Arrange dynamic
sig = sig.*DINAMICA/((2^PRECISION)*GAIN).*1e6; % The data are now converted in microvolt


% for ch = 1 : size(sig,1) - 1
%     sig(ch,:) = remove_int(sig(ch,:), FSamp, 15, 200);    
% end

a(1:12,1) = [12:-1:1]';
a(1:13,2) = [13: 1:25]';
a(1:13,3) = [38:-1:26]';
a(1:13,4) = [39: 1:51]';
a(2:13,5) = [63:-1:52]';

% Longitudinal SD configuration (Added by Ales): Channels 13, 26, 39 and 52 are discarded (they
% correspond to transversal SD configuration - connections among different columns of electrodes)
% a(1:12,1) = -[12:-1:1]';
% a(1:12,2) = [14: 1:25]';
% a(1:12,3) = -[38:-1:27]';
% a(1:12,4) = [40: 1:51]';
% a(2:12,5) = -[63:-1:53]';

% a(1:11,1) = -[11:-1:1]';
% a(1:12,2) = [13: 1:24]';
% a(1:12,3) = -[37:-1:26]';
% a(1:12,4) = [39: 1:50]';
% a(1:12,5) = -[63:-1:52]';



for k1=1:size(a,1)
    for k2=1:size(a,2)
        if (a(k1,k2)==0)            
            sig_out{k1,k2}={};
        else
            sig_out{k1,k2}=sign(a(k1,k2))*sig(abs(a(k1,k2)),:);
        end
    end
end

force = 0*sig(end,:);
%sig_out=sig(1:end-1,:);
