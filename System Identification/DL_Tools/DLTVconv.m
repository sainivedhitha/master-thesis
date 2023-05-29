function [ y ] = DLTVconv(x,fil,nsides)

%$Revision: 217 $
%$Author: dludvig $
%$Date: 2014-08-19 12:04:42 -0500 (Tue, 19 Aug 2014) $

y = zeros(size(x));
if nsides == 1
    nlags = size(fil,1);
else
    nlags = (size(fil,1)+1)/2;
end
ny = length(y);

if nsides == 1
    for i = nlags:ny
        y(i) = x(i-nlags+1:i)'*flipud(fil(:,i));
    end
else
    for i = nlags:ny-nlags
        y(i) = x(i-nlags+1:i+nlags-1)'*flipud(fil(:,i));
    end
end

