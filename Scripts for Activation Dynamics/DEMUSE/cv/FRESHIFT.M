function segt=freshift(seg,teta)

SEG=fft(seg);

f=fftshift([-0.5:1/(length(seg)):0.5-1/(length(seg))]);

SEGt=SEG.*exp(j*2*pi*teta*f);
segt=real(ifft(SEGt));

