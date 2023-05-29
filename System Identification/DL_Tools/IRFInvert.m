function [HC, vafc, FS] = IRFInvert(HI)
R = nldat(randn(10e3,1),'domainincr',.04);
h = irf(cat(2,R,R),'nlags',4,'nsides',2);

% Create low-pass filtered noise for later    
P = parafcn('func',@paratf,'nparams',3);
set(P,'Parameter 1',5);
set(P,'Parameter 2',0.5);
set(P,'Parameter 3',0.15);
set(P,'domainincr',.01);
set(P,'domainmax',4.98);


for i = 1:size(HI,2)
    i
    set(h,'data',-HI(:,i));
    R2 = nlsim(h,R);
    Hc = irf(cat(2,R2,R),'nlags',250);
    HC(:,i) = double(Hc);
    vafc(i) = double(vaf(Hc,cat(2,R2,R)));
    try
        P2 = parafit(P,Hc);
        v = get(P2,'param');
        Vcl(i,:) = cell2mat(get(v,'value'));
        vafr(i) = double(vaf(para2irf(P2),cat(2,R2,R)));       
        vafp(i) = double(vaf(P2,Hc));
    catch
        Vcl(i,:) = [0 0 0];
        vafr(i) = 0;
        vafp(i) = 0;
    end
    F = fresp(cat(2,R,R2));
    FS(:,i) = 20*log10(abs(double(F(:,1))));
    
end

