function [Txy F] = tfestimatecl(X,Y,U,varargin)

[PUY F] = cpsd(U,Y,varargin{:});
[PUX ~] = cpsd(U,X,varargin{:});

Txy = PUY./PUX;