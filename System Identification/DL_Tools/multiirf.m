function [Hident phixx] = multiirf(X,Y,startpoint,endpoint,nlags,nsides,dt)

%$Revision: 33 $
%$Author: dludvig $
%$Date: 2012-05-17 16:26:00 -0500 (Thu, 17 May 2012) $

X = X-mean(mean(X(startpoint:endpoint,:))); %remove offsets
Y = Y-mean(mean(Y(startpoint:endpoint,:)));

if nsides == 1
    m1 = 0;
    m2 = nlags;
  %  Cxx = multicor (X, X, startpoint, endpoint, nlags, 2);
elseif nsides == 2
    m1 = -nlags;
    m2 = nlags;
  %  Cxx = multicor (X, X, startpoint, endpoint, nlags*2, 2);
end


Cyx = multicor (X, Y, startpoint, endpoint, nlags, nsides); % generate multicorrelation matrix

L = numel(X(startpoint:endpoint,:));
hlength = length(Cyx);

phixx = zeros(length(Cyx));
for i = 1:hlength
    Cxx = multicor (X, X, startpoint-m1-i+1, endpoint-m1-i+1, nlags*nsides, 2);
    phixx(:,i) = Cxx(hlength-i+1:end-i+1);
end

% Until here it is explained in the paper

[U,S,V] = svd(phixx);
S = diag(S);
tol = max(size(phixx)')*max(S)*eps;
%tol = 1e-5;
PHIxx_rank = sum(S>tol);


% find Sinv
Sinv = zeros(length(S),1);
Sinv(1:PHIxx_rank) = 1./S(1:PHIxx_rank);
Sinv = diag(Sinv);


out_var = Sinv*((U'*Cyx).^2);
cumul_out_var = cumsum(out_var);
d = (1:hlength)';
yvar = mean(mean(Y(startpoint:endpoint,:).*Y(startpoint:endpoint,:)))-mean(mean(Y(startpoint:endpoint,:)))^2;
MDL = (ones(hlength,1)+log(L)/L*d).*(yvar*ones(hlength,1)-cumul_out_var);
[minMDL, r] = min(MDL);

pseudoinv = V(:,1:r)*Sinv(1:r,1:r)*U(:,1:r)';
Hident = 1/dt*pseudoinv*Cyx;
% Hident = 1/dt*(phixx^-1)*Cyx; 

