% 
close all;
clear all;
clc;

% Calculate valve motion
p = [0.2212, 0.6, 0, 8.566, 1, 0, 0.0433, 0.141, 10.38, 0, 0.8];
%pres = ((p(2))/(p(8)*(9)))^2;
pres = 5;
[t,y] = QuarterWaveSimulation(p,[0.99,0,pres,0,0],[0,1.5],[0,25],'FixedPressure',true);

q = p(8)*p(9)*sqrt(y(3,:)).*y(1,:);
%q(q>2) = 2;

figure;
plot(t,y(1,:))
xlabel('Time - $\tau$','Interpreter','latex');
ylabel('Position - $\tilde{x}$','Interpreter','latex')
xlim([0,6])

figure;
plot(t,q,'r--')
hold on
q(q>2) = 2;
plot(t,q,'b-')
xlabel('Time - $\tau$','Interpreter','latex');
ylabel('Flow rate - $q$','Interpreter','latex')
legend('Required flow rate','Capped flow rate')
xlim([0,6])

figure;
plot(t,y(2,:))
xlabel('Time - $\tau$','Interpreter','latex');
ylabel('Velocity - $\tilde{x}''$','Interpreter','latex')
xlim([0,6])

figure;
plot(t,y(3,:))
xlabel('Time - $\tau$','Interpreter','latex');
ylabel('Tank Pressure - $\tilde{p}$','Interpreter','latex')
xlim([0,6])