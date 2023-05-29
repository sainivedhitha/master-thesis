function [Yp vaft] = multifil(Hident,X,Y,startpoint,endpoint,nsides,dt)

%$Revision: 7 $
%$Author: dludvig $
%$Date: 2012-02-20 10:25:42 -0600 (Mon, 20 Feb 2012) $

y_temp =[];
Yp =[];
for i =1:size(Y,2) %1:nr
    if nsides == 1
        y_temp(:,1) = dt*DLconv(X(:,i),Hident,1);
    elseif nsides == 2
        y_temp(:,1) = dt*DLconv(X(:,i),Hident,2);
    end
    Yp(:,i) = y_temp(startpoint:endpoint,1); %estimated output
end

Yr = Y(startpoint:endpoint,:)-Yp; %residuals

varYr = mean(mean(Yr.^2))-(mean(mean(Yr))).^2;
varY = mean(mean(Y(startpoint:endpoint,:).^2))-(mean(mean(Y(startpoint:endpoint,:)))).^2;

vaft = 100*(1-varYr/varY);