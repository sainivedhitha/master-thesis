function [A Hident Yp vaft] = multinl2(X,Y,startpoint,endpoint,nlags,nsides,norder,dt)

%$Revision: 7 $
%$Author: dludvig $
%$Date: 2012-02-20 10:25:42 -0600 (Mon, 20 Feb 2012) $

Hident = multiirf(X,Y,startpoint,endpoint,nlags,nsides,dt);
% Hident = zeros(nlags+1,1);
% Hident(1) = 1;


Hidentnew = 0;
Anew = 0;
Ypnew = 0;
vaftnew = 0;
delta_vaf = 1;
vaf_old = 0;
while delta_vaf > .05
    PX = [];
    
    if nsides == 1
        j1 = 0;
        j2 = nlags;
    elseif nsides == 2
        j1 = -nlags;
        j2 = nlags;
    end
    
    b = max(max(X(startpoint-j2:endpoint-j1)));
    a = min(min(X(startpoint-j2:endpoint-j1)));
    
    Xnorm = (2*X-(b+a))/(b-a);
    
    for r = 1:size(X,2)
        for ti = startpoint:endpoint
            x_temp = DLtcheb(Xnorm(ti-j2:ti-j1,r),norder);
            for n=1:norder+1
                px(ti-startpoint+1,n) = sum(x_temp(:,n).*flipud(Hident));
            end
        end
        PX = [PX; px];
    end
    
    y_temp = reshape(Y(startpoint:endpoint,:),numel(Y(startpoint:endpoint,:)),1);
    pseudoinv = pinv2(PX);
    A = pseudoinv*y_temp;
    
    for i = 1:size(X,2)
        X2(:,i) = chebval(A,X(:,i),a,b);
    end
    
    Hident = multiirf(X2,Y,startpoint,endpoint,nlags,nsides,dt);
    [Yp vaft] = multifil(Hident,X2,Y,startpoint,endpoint,nsides,dt);
    delta_vaf = vaft-vaf_old;
    vaf_old = vaft;
    gain=sum(Hident)*dt;
    Hident = Hident/gain;
    A = A*gain;
    disp(['Hammerstein VAF = ' num2str(vaft) '%']);
    if delta_vaf > 0
        Hidentnew = Hident;
        Anew = A;
        Ypnew = Yp;
        vaftnew = vaft;
    end
end
Hident = Hidentnew;
A = Anew;
Yp = Ypnew;
vaft = vaftnew;
gain=sum(Hident)*dt;
Hident = Hident/gain;
A = A*gain;
end
