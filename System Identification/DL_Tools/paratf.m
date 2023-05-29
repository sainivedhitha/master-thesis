function y = paratf (t, params)
% Function for non-linear regression to any transfer function

%$Revision: 7 $
%$Author: dludvig $
%$Date: 2012-02-20 10:25:42 -0600 (Mon, 20 Feb 2012) $

% Define your transfer function or use global sys

%sys = tf(params(1)*params(3)^2,[1 2*params(2)*params(3) params(3)^2]);
sys = tf(1,[params(3) params(2) params(1)]);
% sys = tf(1,[params(3) params(2) params(1)])*tf(1,[1 params(4)]);

% Add other params and options
model_type = 'irf'; %irf, fr, step
% delay = 7;
%sys.inputd = abs(params(4));
nsmo = 1;

switch model_type
    case 'irf'
        Q = impulse(sys,t);
%         Q2 = [zeros(delay,1); Q];
%  y = smo(Q,nsmo);
y=Q;
 
    case 'fr'
        y = freqresp(sys,t);
    case 'step'
        Q = step(sys,t);
%         Q2 = [zeros(delay,1); Q];
        y = smo(Q,nsmo);
end