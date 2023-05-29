function [Hi, Ar, Hr, Yt, Yi, Yr, vaft, Hiv, Hrv] = multipc (X,Y,V,startpoint,endpoint,reflex_lags,norder,dt,dr)

%$Revision: 7 $
%$Author: dludvig $
%$Date: 2012-02-20 10:25:42 -0600 (Mon, 20 Feb 2012) $

[Xd Vd Yd] = deal(zeros(size(X,1)/dr,size(X,2)));
% V = zeros(size(X));

if dr ~= 0 && dr ~= 1
    for i = 1:size(X,2)
        V(:,i) = double(ddt(nldat(X(:,i),'domainincr',.001)));
        Xd(:,i) = decimate(X(:,i),dr);
        Vd(:,i) = decimate(V(:,i),dr);
        Yd(:,i) = decimate(Y(:,i),dr);
    end
else
    Xd = X;
    Yd = Y;
    Vd = V;
end

Yi = zeros(size(Yd));
Yr = zeros(size(Yd));
Yrr = Yd;

Ytnew = 0;
Hinew = 0;
Arnew =0;
Hrnew = 0;
vaftnew = 0;
delta_vaf = 1;
vaf_old = 0;
while delta_vaf > .05
    [Hi pinv] = multiirf(Xd,Yrr,startpoint,endpoint,.04/(dt*dr),2,dt*dr);
    Yp = multifil(Hi,Xd,Yd,startpoint,endpoint,2,dt*dr);
    Yi(startpoint:endpoint,:) = Yp;
    Yir = Yd-Yi;
    Hiv = diag(pinv)*mean(var(Yrr(startpoint:endpoint,:)-Yp));
    
%     [Ar Hr Yp ~] = multinl2(Vd,Yir,startpoint,endpoint,reflex_lags,1,norder,dt*dr);
Ar = 0;
    [Hr pinv] = multiirf(Vd,Yir,startpoint,endpoint,reflex_lags,1,dt*dr);
    Yp = multifil(Hr,Vd,Yir,startpoint,endpoint,1,dt*dr);
    Yr(startpoint:endpoint,:) = Yp;
    Yrr = Yd-Yr;
    Yt = Yi+Yr;
    
    Hrv = diag(pinv)*mean(var(Yir(startpoint:endpoint,:)-Yp));
    
    Ytr = Yd-Yi-Yr;
    varYtr = mean(mean(Ytr(startpoint:endpoint,:).^2))-(mean(mean(Ytr(startpoint:endpoint,:)))).^2;
    varY = mean(mean(Yd(startpoint:endpoint,:).^2))-(mean(mean(Yd(startpoint:endpoint,:)))).^2;
    
    vaft = 100*(1-varYtr/varY);
    delta_vaf = vaft-vaf_old;
    vaf_old = vaft;
    disp(['Parralel-Cascade VAF = ' num2str(vaft) '%']);
        if delta_vaf > 0
            Hinew = Hi;
            Arnew =Ar;
            Hrnew = Hr; 
            vaftnew = vaft;
            Ytnew = Yt;
        end 
end

Yt = Ytnew;
Hi = Hinew;
Ar = Arnew;
Hr = Hrnew;
vaft = vaftnew;

end
