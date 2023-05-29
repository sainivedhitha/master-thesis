function mvar = multivar(x)

%$Revision: 7 $
%$Author: dludvig $
%$Date: 2012-02-20 10:25:42 -0600 (Mon, 20 Feb 2012) $

 mvar = mean(mean(x.^2))-(mean(mean(x))^2);


end

