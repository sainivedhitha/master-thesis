function C = multicor (X, Y, startpoint, endpoint, nlags, nsides)
% Generate multicorrelation matrix C=phixy. 

%$Revision: 7 $
%$Author: dludvig $
%$Date: 2012-02-20 10:25:42 -0600 (Mon, 20 Feb 2012) $

if nsides == 1
    y = Y(startpoint:endpoint,:);
    for i = 1:nlags+1
        x = X(startpoint+1-i:endpoint+1-i,:);
        phixy = mean(mean(x.*y)); %average over time and repetitions
        C(i,1) = phixy;
    end
elseif nsides == 2
     y = Y(startpoint:endpoint,:);
     for i = 1:nlags*2+1
        x = X(startpoint+nlags+1-i:endpoint+nlags+1-i,:);
        phixy = mean(mean(x.*y));
        C(i,1) = phixy;
     end
end
end

