function [sig, vpp, epocs] = adjst2(sig1, sig2, sig3, sig4,FSamp);

sig = [sig1(1 : 15,1 : end); sig2(1 : 15,1 : end); sig3(1 : 15,1 : end); sig4(1 : 15,1 : end)];
epocs = floor(size(sig1,2)/FSamp);
[Row Col] = size(sig);
temp = 0;
vpp = 0;
for column = 1 : Col
   for row = 1 : Row
      temp = max(max(abs(sig(row,column))))/1.5;
      if temp > vpp
         vpp = temp;
      end
   end
end