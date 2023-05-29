function teta=localac(sig1,sig2,dint,fsamp)

% Search limits for localac

CVmin=2;
CVmax=60;  % prej 7
TETAmin=floor(dint/CVmax*fsamp);
TETAmax=ceil(dint/CVmin*fsamp);

for i = TETAmin : TETAmax
    corrloc(i-TETAmin+1)=sum(sig1(1:length(sig1)-i).*sig2(i+1:length(sig2)));
end;

delay=[TETAmin:TETAmax];

[a b]=max(corrloc);

if b > 1 & b < length(delay),
    
    [P S]=polyfit(delay([b-1 b b+1]),corrloc([b-1 b b+1]),2);
    
    teta=-P(2)/(2*P(1));
    
else,
    
    teta=b;
    
end;


