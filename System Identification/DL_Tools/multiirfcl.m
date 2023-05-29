function [Hident pseudoinv] = multiirfcl(X,Y,U,startpoint,endpoint,nlags,nsides,dt)
% Multi-segment closed-loop IRF
% X: input
% Y: Output
% U: Intsrumental Variable
% startpoint: start of segment must be at least 2x nlags from start
% endpoint: end of segment must be at least 2x nlags from end
% nlags: number of lags (1 side) in IRF
% nsides: number of sides
% dt: sampling increment
% X
% return;

%$Revision: 34 $
%$Author: dludvig $
%$Date: 2012-05-17 16:33:07 -0500 (Thu, 17 May 2012) $

X = X-mean(mean(X(startpoint:endpoint,:)));
Y = Y-mean(mean(Y(startpoint:endpoint,:)));
U = U-mean(mean(U(startpoint:endpoint,:)));

if nsides == 1
    m1 = 0;
    m2 = nlags;
  %  Cxx = multicor (X, X, startpoint, endpoint, nlags, 2);
elseif nsides == 2
    m1 = -nlags;
    m2 = nlags;
  %  Cxx = multicor (X, X, startpoint, endpoint, nlags*2, 2);
end
Cyu = multicor (U, Y, startpoint, endpoint, nlags, nsides);

L = numel(X(startpoint:endpoint,:));
hlength = length(Cyu);

phixx = zeros(length(Cyu));
for i = 1:hlength
    Cxu = multicor (U, X, startpoint-m1-i+1, endpoint-m1-i+1, nlags*nsides, 2);
    phixu(:,i) = Cxu(hlength-i+1:end-i+1);
end

[Um,S,V] = svd(phixu);
S = diag(S);
tol = max(size(phixu)')*max(S)*eps;
PHIxu_rank = sum(S>tol);

% find Sinv
Sinv = zeros(length(S),1);
Sinv(1:PHIxu_rank) = 1./S(1:PHIxu_rank);
Sinv = diag(Sinv);


out_var = Sinv*((Um'*Cyu).^2);
cumul_out_var = cumsum(out_var);
d = (1:hlength)';
yvar = mean(mean(Y(startpoint:endpoint,:).*Y(startpoint:endpoint,:)))-mean(mean(Y(startpoint:endpoint,:)))^2;
MDL = (ones(hlength,1)+log(L)/L*d).*(yvar*ones(hlength,1)-cumul_out_var);
[minMDL, r] = min(MDL);

pseudoinv = V(:,1:r)*Sinv(1:r,1:r)*Um(:,1:r)';
Hident = 1/dt*pseudoinv*Cyu;
% Hident = 1/dt*(phixx^-1)*Cyx; 

