function [cv,teta] = mle3(Segna,start,dint,fsamp)

cpter=1;
t=start;
teta=10;
trial=0; 
num_sig=size(Segna,1);

while (abs(teta-t)>=5e-5) & trial < 30
    trial=trial+1;
    teta=t;     
    de1=0;
    de2=0;
    for i=1:num_sig,
        
        [de1t, de2t]=derivbeam(Segna,i,teta);
        
        de1=de1+de1t+eps;
        
        de2=de2+de2t+eps;
        
    end;
    
    % Newton's criteria
    if (de2>0) 
        u=-de1/de2;
        if (abs(u)>0.5)
            u=-0.5*abs(de1)/de1;
        end
    else
        u=-0.5*abs(de1)/de1;
    end
    %err1(cpter,1)=de1;
    %err2(cpter,1)=de2;
    %u_v(cpter,1)=u;
    %cpter=cpter+1;
    
    %u=-de1/de2;
    
    t=teta+u;	  % result
end

cv=dint/(teta/fsamp);