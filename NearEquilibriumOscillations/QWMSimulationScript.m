% Reset
close all;
clear all;
clc;

model = 'QWM';
valve = 'pilot';

% Parameters to change
cont = false; % set to true to simulate after initial collision
q = 0.6;
L = 20; % 0.3, 0.35, 0.5
x_range = [0,1.1];
t_range = [0,20];
%t_range = t_range(1):(5e-3):t_range(2);
%t_range = linspace(t_range(1),t_range(2),4000);

% Set dimensionless parameter vector
switch valve
    case 'pilot'
        p = CalculateParameters(q,L);
    case 'spring'
        p = [L*201.3/890, q, (1.78e-3)/L, 3.186, 1, 0, 0.0141, 0.0538, 2.352, 0.00064]; % Dimensionless spring PRV 1E2 parameters
end
%p(7) = 0;
p(3) = 0; p(10) = 0; %Lamdba must equal 0, but phi = 0 has little effect

pres = ((p(2))/(p(8)*p(9)))^2;
init = [0.99,0,pres,0,0]; % need 10000 to see some oscillations

% Calculate QWM
[t_control, y_control] = QuarterWaveSimulation(p,init,x_range,t_range,'Initial',cont);
[t,y] = QuarterWaveSimulation(p,init,x_range,t_range,model,cont);
%[t2,y2] = QuarterWaveSimulation(p,init,x_range,t_range,model);

% Plot graphs
figure;
plot(t_control,y_control(1,:),'r--')
hold on
plot(t,y(1,:),'b-')
%plot(t2,y2(1,:))
xlabel('Time'); ylabel('Valve position')
legend('Initial','Quarter Wave Model','location','northeast')

figure;
plot(t_control,y_control(2,:),'r--')
hold on
plot(t,y(2,:),'b-')
%plot(t2,y2(2,:))
xlabel('Time'); ylabel('Valve velocity')
legend('Initial','Quarter Wave Model','location','northeast')

figure;
plot(t_control,y_control(3,:),'r--')
hold on
plot(t,y(3,:),'b-')
%plot(t2,y2(3,:))
xlabel('Time'); ylabel('Tank pressure')
legend('Initial','Quarter Wave Model','location','northeast')

figure;
plot(t_control,y_control(4,:),'r--')
hold on
plot(t,y(4,:),'b-')
%plot(t2,y2(4,:))
xlabel('Time'); ylabel('Pressure fluctuation')
legend('Initial','Quarter Wave Model','location','northwest')

figure;
plot(t_control,y_control(5,:),'r--')
hold on
plot(t,y(5,:),'b-')
%plot(t2,y2(5,:))
xlabel('Time'); ylabel('Velocity fluctuation')
legend('Initial','Quarter Wave Model','location','northwest')