function [de1, de2]=derivbeam(Segna,posref,teta)

%Calculation of first and second devivative of the mean square error 
%for the beamforming technique in the case of first signal as reference

num_sig=size(Segna,1);

M=num_sig-1;

N=size(Segna,2);

pos=find([1:M+1]~=posref);

Htemp=[Segna(posref,:);Segna(pos,:)];

Segna=Htemp;

position(1)=posref;
for i=2:posref
   
   position(i)=i-(posref+1);
   
end;

for i=1:num_sig-posref
   
   position(i+posref)=i;
   
end;

position=position(2:num_sig);

k=[1:N/2];

terminede1=zeros(1,length(k));

terminede2=zeros(1,length(k));

terminede12=zeros(1,length(k));

terminede22=zeros(1,length(k));

for i=1:num_sig
   SEGNAfft(i,:)=fft(Segna(i,:));
end;

%Calcolo della derivata prima

for i=1:M
   for u=i+1:M
   
         terminede1=terminede1-imag(SEGNAfft(i+1,k+1).*conj(SEGNAfft(u+1,k+1)).*exp(j*2*pi*k*(position(i)-position(u))*teta/N).*2*pi.*k.*(position(i)-position(u))./N);
         
   end;
end;
terminede1=terminede1*2/(M^2);

for i=1:M
   
   terminede2=terminede2+SEGNAfft(i+1,k+1).*exp(j*2*pi*k*position(i)*teta/N).*2*pi.*k.*position(i)./N;

end;

terminede2=2*imag(conj(SEGNAfft(1,k+1)).*terminede2)/M;

de1=2/N*sum(terminede1+terminede2);

%Calcolo della derivata seconda

for i=1:M
   for u=i+1:M
   
         terminede12=terminede12-real(SEGNAfft(i+1,k+1).*conj(SEGNAfft(u+1,k+1)).*exp(j*2*pi*k*(position(i)-position(u))*teta/N).*((2*pi*k*(position(i)-position(u))/N).^2));
         
   end;
end;
terminede12=terminede12*2/(M^2);

for i=1:M
   
   terminede22=terminede22+SEGNAfft(i+1,k+1).*exp(j*2*pi*k*position(i)*teta/N).*((2*pi*k*position(i)/N).^2);

end;

terminede22=2*real(conj(SEGNAfft(1,k+1)).*terminede22)/M;

de2=2/N*sum(terminede12+terminede22);


