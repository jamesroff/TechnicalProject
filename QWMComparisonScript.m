% Reset
close all;
%clear all;
%clc;

valve = 'pilot';

% Parameters to change
cont = false; % set to true to simulate after initial collision

q = 0.6;
%%% Very weird behaviour with L = 9 - L = 10 - L = 11 for 'spring'
L = 3; % 0.3, 0.35, 0.5
x_range = [0,1.1];
t_range = [0,20];
%t_range = t_range(1):(5e-3):t_range(2);
%t_range = linspace(t_range(1),t_range(2),4000);

% Set dimensionless parameter vector
switch valve
    case 'pilot'
        p = CalculateParameters(q,L);
    case 'spring'
        p = [L*201.3/890, q, (1.78e-3)/L, 3.186, 1, 0, 0.0141, 0.0538, 2.352, 0.00064, 0.8]; % Dimensionless spring PRV 1E2 parameters
end
%p(7) = 0; % beta = 0;
%p(3) = 0; p(10) = 0; % Lambda = 0; phi = 0;

pres = ((p(2))/(p(8)*p(9)))^2;
init = [0.8,0,pres,0,0];

% Calculate QWM
[t1, y1] = QuarterWaveSimulation(p,init,x_range,t_range,'FullQWM',cont);
[t2, y2] = QuarterWaveSimulation(p,init,x_range,t_range,'QWM',cont);

% Plot graphs
figure;
plot(t1,y1(1,:),'r--')
hold on
plot(t2,y2(1,:),'b-')
%plot(t2,y2(1,:))
xlabel('Time'); ylabel('Valve position')
legend('Initial','Quarter Wave Model','location','northeast')

figure;
plot(t1,y1(2,:),'r--')
hold on
plot(t2,y2(2,:),'b-')
%plot(t2,y2(2,:))
xlabel('Time'); ylabel('Valve velocity')
legend('Initial','Quarter Wave Model','location','northeast')

figure;
plot(t1,y1(3,:),'r--')
hold on
plot(t2,y2(3,:),'b-')
%plot(t2,y2(3,:))
xlabel('Time'); ylabel('Tank pressure')
legend('Initial','Quarter Wave Model','location','northeast')

figure;
plot(t1,y1(4,:),'r--')
hold on
plot(t2,y2(4,:),'b-')
%plot(t2,y2(4,:))
xlabel('Time'); ylabel('Pressure fluctuation')
legend('Initial','Quarter Wave Model','location','northwest')

figure;
plot(t1,y1(5,:),'r--')
hold on
plot(t2,y2(5,:),'b-')
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