function [p] = CalculateParameters(q,L)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Define physical parameters
m_v = 0.4392;
c_v = 0;
C_d = 0.32;
D = 0.0266;
A_p = (pi*(D^2))/4;
A_v = 1.5*A_p;
m_cap = 3.8;
rho = 998.2;
a = 890;
V = 10.6;
r = 0.8;
lambda = 0.02;

zeta = sqrt(2)*pi*D*C_d;

% Calculate reference parameters
x_ref = (1/C_d)*sqrt((A_v-A_p)/(4*pi));
p_ref = 1e5;
w_ref = C_d * sqrt(4*pi) * sqrt((p_ref * x_ref)/(m_v));

% Calculate dimensionless parameters
gamma = (L*w_ref)/(a);
%Lambda = x_ref / L;
Lambda = 0;
alpha = (w_ref*a*rho*x_ref)/(p_ref);
delta = (zeta*x_ref)/(A_p);
kappa = (c_v)/(m_v*w_ref);
beta = ((a^2)/(V))*((m_cap)/(p_ref*w_ref));
%beta = 0;
mu = (rho*A_p*w_ref*x_ref)/(m_cap);
sigma = (zeta*sqrt(p_ref))/(A_p*sqrt(rho)*w_ref);
%phi = lambda*((x_ref)/(2*D));
phi = 0;

%omega = 

p = [gamma, q, Lambda, alpha, delta, kappa, beta, mu, sigma, phi, r];

end