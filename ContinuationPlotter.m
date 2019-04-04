%% 
close all;
clear all;

%% Load continuation data calculated in MatCont
load H_H(1).mat;
p_right = x(5,:);
gamma_right = x(6,:);
clearvars -except p_right gamma_right

load H_H(2).mat;
p_left = x(5,:);
gamma_left = x(6,:);
clearvars -except p_right gamma_right p_left gamma_left

p_cont = [p_left(end:-1:2), p_right];
gamma_cont = [gamma_left(end:-1:2), gamma_right];
clearvars -except p_cont gamma_cont
clc

%% Calculate analytic stability boundary
kappa = 0; alpha = 8.5658; delta = 1; sigma = 10.3808;

gamma = linspace(0,20,2001); gamma(1) = [];

pres = ((sqrt(2)/2)*alpha*(pi./(2*gamma)).^3)./(delta^2 - alpha*(delta^2)*(pi./(2*gamma)));
gamma = gamma(pres>0); pres = pres(pres>0);

%% Plot two dimensional bifurcation diagram
figure;
plot(p_cont,gamma_cont)
hold on
plot(pres,gamma,'r--')
plot(0.0703,14.9501,'rx')
plot(0.1686,1.4745,'b*')
ylabel('\gamma'); xlabel('p_0')
xlim([0,1]); ylim([0,20])
grid on