% Reset
close all;
clear all;
clc;

% Set dimensionless parameter vector
q = 0.6;
L = 0.35; % 0.3
p = CalculateParameters(q,L);

% Calculate QWM
pres = ((p(2))/(p(8)*p(9)))^2;
init = [0.9,0,1,0,0];
x_range = [0,1.1];
t_range = [0,10];
%
%t_range = linspace(0,20,10000);
%
model = 'QWM';
[t,y] = QuarterWaveSimulation(p,init,x_range,t_range,model,0);
%[t2,y2] = QuarterWaveSimulation(p,init,x_range,t_range,model);

% Plot graphs
figure;
plot(t,y(1,:))
%plot(t2,y2(1,:))
xlabel('Time'); ylabel('Valve position')

figure;
plot(t,y(2,:))
%plot(t2,y2(2,:))
xlabel('Time'); ylabel('Valve velocity')

figure;
plot(t,y(3,:))
%plot(t2,y2(3,:))
xlabel('Time'); ylabel('Tank pressure')

figure;
plot(t,y(4,:))
%plot(t2,y2(4,:))
xlabel('Time'); ylabel('Pressure fluctuation')

figure;
plot(t,y(5,:))
%plot(t2,y2(5,:))
xlabel('Time'); ylabel('Velocity fluctuation')