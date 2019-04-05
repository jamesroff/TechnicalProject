% Reset
close all;
%clear all;
%clc;

model = 'QWM';
valve = 'pilot';

% Parameters to change
cont = false; % set to true to simulate after initial collision

q = 0.6;
%%% Very weird behaviour with L = 9 - L = 10 - L = 11 for 'spring'
L = 20; % 0.3, 0.35, 0.5
x_range = [0,1.1];
t_range = [0,100];
%t_range = t_range(1):(5e-3):t_range(2);
%t_range = linspace(t_range(1),t_range(2),4000);

% Set dimensionless parameter vector
switch valve
    case 'pilot'
        p = CalculateParameters(q,L);
    case 'spring'
        p = [L*201.3/890, q, (1.78e-3)/L, 3.186, 1, 0, 0.0141, 0.0538, 2.352, 0.00064, 0.8]; % Dimensionless spring PRV 1E2 parameters
end
p(7) = 0; % beta = 0;
p(3) = 0; p(10) = 0; % Lambda = 0; phi = 0;

a = 0.9;
p(1) = (pi/2)*p(4) / a;

CalculateFrequency(p(1));

%pres = ((p(2))/(p(8)*p(9)))^2;
%pres = ((sqrt(2)*pi/4)*(p(4)/p(1))*(pi/(2*p(1)))^2)/(p(5)^2 * 2 * (1-a)); %+ 0.07;
pres = 1;
%pres = 10 * pres;
%pres = pres / 10;
init = [1,0,pres,-0.001,0]; % need 10000 to see some oscillations
%init = [0.99,0,pres,0,0];
%init = [0.99,0,1,0,0];
%init = [1,0,1,0.01,0];

% Calculate QWM
[t_control, y_control] = QuarterWaveSimulation(p,init,x_range,t_range,'Initial',true);
[t,y] = QuarterWaveSimulation(p,init,x_range,t_range,model,cont);
%[t2,y2] = QuarterWaveSimulation(p,init,x_range,t_range,model);
if isempty(t) && isempty(y)
    fprintf("Error during Quarter Wave Simulation\n");
    return
end

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


% %% FFT
% 
% T = t(3):5e-3:t(end);
% Y = interp1(t(3:end),y(:,3:end)',T); Y = Y';
% % figure;
% % plot(T,Y(1,:))
% 
% FY = fft(Y'); FY = FY';
% 
% % FY = fft(y'); FY = FY';
% 
% P2 = abs(FY/size(FY,2));
% P1 = P2(:,1:floor(size(FY,2)/2)+1);
% P1 = 2*P1;
% 
% Fs = 1/(5e-3); % sampling frequency
% f = Fs*(0:floor(size(FY,2)/2))/size(FY,2);
% 
% figure;
% semilogy(f,P1(1,:))