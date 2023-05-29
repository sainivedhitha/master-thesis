function [y] = DLconv(x,fil,nsides)

%$Revision: 7 $
%$Author: dludvig $
%$Date: 2012-02-20 10:25:42 -0600 (Mon, 20 Feb 2012) $

X = zeros(length(x),length(fil));

if nsides == 1

for i = 1:length(fil)
    X(:,i) = [zeros(i-1,1); x(1:end+1-i)]; 
end

elseif nsides == 2
    hlength = (length(fil)-1)/2;
    for i = 1:hlength
        X(:,hlength-i+1) = [x(i+1:end); zeros(i,1)];
    end
    for i = 1:hlength+1
        X(:,hlength+i) = [zeros(i-1,1); x(1:end+1-i)];
    end
end

y = X*fil;

end

