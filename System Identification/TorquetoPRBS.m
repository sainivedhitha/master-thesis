% TorquetoPRBS transforms the torque perturbation into a PRBS position
% signal. 

clear all
close all
clc

% Load the OutputData from Achilles
load('D:\SaiArunajatesan\Thesis\Data\subject06\Trial3C1R8');

Torque = OutputData.Drivefile;
fs = 2048;                                  % Sampling frequency [Hz]
dt = 1/fs;                                  % Sampling time [s]
t = 0:1/fs:length(Torque)/fs-1/fs;            % Time vector [s]
Tmax = max(t);                              % Total duration of the pert. [s]

% Virtual parameters of the admittance controller of the Achilles
% Virtual parameters of the admittance controller of the Achilles
I = 0.1;                                      % Inertia [kgm^2]
D = 2.5;                                    % Damping [Nm/rad/s]
K = 4500;                                     % Stiffness [Nm/rad]


%%
B = [I D K];
A = 1;

% In order to multiply the inverse of the admittance controller with the
% low pass filter, it is necessary to transform the transfer function (in
% the Laplace domain) to the Z-transform. This is done with the command c2d
H = tf(B,A);
H_inv = H^-1;
H_inv_z = c2d(H_inv, dt, 'tustin');

% maxflat command allows to design filters with the desired number of poles
% and zeros. In this case, the filter should have two more poles than zeros
% (to have a proper TF).
fc = 40;
[D,C] = maxflat(0,2,2*fc*dt);

H_filter_z = tf(D,C,dt);
H_final_z = H_filter_z*H_inv_z;

% Just in case the continuous coefficients want to be checked. 
H_filter_c = d2c(H_filter_z);
H_final_c = H_filter_c*H_inv;

% Obtaining the parameters of the filter (discrete domain because 'filter'
% command works with coefficients in this domain)
B = H_final_z.Numerator;
B = B{1};
A = H_final_z.Denominator;
A = A{1};

PRBS_trick = filter(B,A,Torque);
% figure, plot(PRBS_trick)

% Save PRBS data
dir = strcat(strcat('D:\SaiArunajatesan\Thesis\Data\subject06\PRBS_6-8','.mat'));
save(dir, 'PRBS_trick');