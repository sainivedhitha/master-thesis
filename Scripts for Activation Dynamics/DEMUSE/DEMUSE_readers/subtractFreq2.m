function ynew=subtractFreq2(y,Ind,r,draw)
% subtracts the content at the frequencies [Ind-r:Ind+r] from the signal y. With paramter draw>0 it also plots
% the frequency content of original and filtered signal.
% OUTPUT
%   ynew - filtered signal

fy=fft(y);
lenFy=length(fy);
fd=zeros(1,lenFy);  
tInd=[];
% for k=-floor(r/2):floor(r/2)
%     tInd=[tInd Ind+k];
% end

% mFreq=median(fy);
% sFreq=std(fy);
% t2Ind=find(abs(fy) > mFreq+15*sFreq);
% for k=-floor(r/2):floor(r/2)
%     tInd=[tInd t2Ind+k];
% end

WinLen=1000;
for k1=1:WinLen:length(fy)-WinLen;
    mFreq=median(abs(fy(k1+1:k1+WinLen)));
    sFreq=std(abs(fy(k1+1:k1+WinLen)));
    t2Ind=find(abs(fy(k1+1:k1+WinLen)) > mFreq+5*sFreq);
    t2Ind=t2Ind+k1;
    for k=-floor(r/2):floor(r/2)
        tInd=[tInd t2Ind+k];
    end
end


tInd=round(tInd(find(tInd>0 & tInd<=lenFy/2+1)));
for k=tInd
    fd(k)=fy(k);    
end

correct=lenFy-floor(lenFy/2)*2;
fd(lenFy:-1:ceil(lenFy/2)+1)=conj(fd(2:1:ceil(lenFy/2)+1-correct));

ynew=y-ifft(fd(1:lenFy));

if draw>0
fy=fftshift(fy);
fd=fftshift(fd);     
    figure(draw); 
    subplot(2,1,1); plot(-ceil(lenFy/2)+1:ceil(lenFy/2)-correct,real(fy)); hold on; plot(-ceil(lenFy/2)+1:ceil(lenFy/2)-correct,real(fy-fd),'r'); hold off; legend('original signal','filtered signal'); title('real');     
    subplot(2,1,2); plot(-ceil(lenFy/2)+1:ceil(lenFy/2)-correct,imag(fy)); hold on; plot(-ceil(lenFy/2)+1:ceil(lenFy/2)-correct,imag(fy-fd),'r'); hold off; legend('original signal','filtered signal'); title('imag'); 
    %pause;
end