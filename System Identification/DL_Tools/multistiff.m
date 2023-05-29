function [Hi, Yi, vaft, Hiv] = multistiff(X,Y,startpoint,endpoint,dt,dr,nlags)

%Compute 2 sided stiffness IRF using multiple sgements
%Usage [Hi, Yi, vaft] = multistiff(X,Y,startpoint,endpoint,dt,dr,nlags)
%Inputs:
%         X: input (Usually Position)
%         Y: output (usually Torque)
%         startpoint: timepoint where segment of interest begings (must be 3x nlags from start)
%         endpoint: timepoint where segment of interest ends (must be 3x nlags from end)
%         dt: sampling increment
%         dr: decimation ratio (0 or 1 if no decimation desired)
%         nlags: number of lags on each side (IRF will be computed from -nlags:nlags)
% Outputs:
%         Hi: estimated IRF
%         Yi: predicted output based on IRF
%         vaft: variance accounted for by IRF

%$Revision: 7 $
%$Author: dludvig $
%$Date: 2012-02-20 10:25:42 -0600 (Mon, 20 Feb 2012) $

%% predifine sizes
[Xd Yd] = deal(zeros(size(X,1)/dr,size(X,2))); 

if dr ~= 0 && dr ~= 1
    for i = 1:size(X,2)
        Xd(:,i) = decimate(X(:,i),dr);
        Yd(:,i) = decimate(Y(:,i),dr);
    end
else
    Xd = X;
    Yd = Y;
end

Yi = zeros(size(Yd));
%%
[Hi phixx] = multiirf(Xd,Yd,startpoint,endpoint,nlags,2,dt*dr); %find IRF
%%
Yp = multifil(Hi,Xd,Yd,startpoint,endpoint,2,dt*dr); %simulate IRF
Yi(startpoint:endpoint,:) = Yp;
Yir = Yd-Yi; %residuals

varYir = mean(mean(Yir(startpoint:endpoint,:).^2))-(mean(mean(Yir(startpoint:endpoint,:)))).^2;
varY = mean(mean(Yd(startpoint:endpoint,:).^2))-(mean(mean(Yd(startpoint:endpoint,:)))).^2;

vaft = 100*(1-varYir/varY);

Hiv = diag(phixx)*varYir;

end