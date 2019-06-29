%% Reset
close all;
clear all;
clc;

% Set dimensionless parameter vector
q = 0.6;
L = 0.6; % 0.3, 0.35, 0.5
p = CalculateParameters(q,L);
%p = [L*201.3/890, q, (1.78e-3)/L, 3.186, 1, 0, 0.0141, 0.0538, 2.352, 0.00064]; % Dimensionless spring PRV 1E2 parameters
p(3) = 0; % Set Lambda to 0
p(10) = 0; % Set phi to 0

% Initialise equilibrium pressure and 
pres = ((p(2))/(p(8)*p(9)))^2;
t_range = [0,1];

%% Ideal behaviour
init = [1,0,pres,0.1,0]; % need 10000 to see some oscillations
[t1, y1] = QuarterWaveSimulation(p,init,[0,1.1],t_range,'QWM',false);

% CHANGE BACK TO INCLUDING PROBLEM TERM

%% Actual behaviour
init = [1,0,pres,0.1,0]; % need 10000 to see some oscillations
[t2, y2] = QuarterWaveSimulation(p,init,[0,1.1],t_range,'QWM',false);

%% High pressure
init = [1,0,10000,0.1,0]; % need 10000 to see some oscillations
[t3, y3] = QuarterWaveSimulation(p,init,[0,1.1],t_range,'QWM',false);

%% Plot B
figure;
plot(t1,y1(4,:),'g-.')
hold on
plot(t2,y2(4,:),'r--')
plot(t3,y3(4,:),'b-')
xlabel('Time'); ylabel('Pressure')
legend('Ideal Wave','Real Wave','High Pressure','location','southeast')

%% Plot C
figure;
plot(t1,y1(5,:),'g-.')
hold on
plot(t2,y2(5,:),'r--')
plot(t3,y3(5,:),'b-')
xlabel('Time'); ylabel('Velocity')
legend('Ideal Wave','Real Wave','High Pressure','location','northeast')