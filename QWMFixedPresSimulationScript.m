% Reset
close all;
%clear all;
%clc;

model = 'QWM';
valve = 'pilot';
graph = 'sep'; % 'sep' or 'all'

% Parameters to change
cont = true; % set to true to simulate after initial collision

q = 0.6;
%%% Very weird behaviour with L = 9 - L = 10 - L = 11 for 'spring'
L = 20; % 0.3, 0.35, 0.5
x_range = [0,1.1];
t_range = [0,600]; % 14.1
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

%p(1) = 4; % 2 and 4

CalculateFrequency(p(1));

%pres = ((p(2))/(p(8)*p(9)))^2;
pres = ((sqrt(2)*pi/4)*(p(4)/p(1))*(pi/(2*p(1)))^2)/(p(5)^2 * 2 * (1-a)); %+ 0.07;
%pres = 0.1;
%pres = 1;
init = [1,0,pres,-0.001,0]; % need 10000 to see some oscillations
%init = [0.99,0,pres,0,0];

% Calculate QWM
%[t_control, y_control] = QuarterWaveSimulation(p,init,x_range,t_range,'Initial',true);
[t,y] = QuarterWaveSimulation(p,init,x_range,t_range,model,cont);
%[t2,y2] = QuarterWaveSimulation(p,init,x_range,t_range,model);
if isempty(t) && isempty(y)
    fprintf("Error during Quarter Wave Simulation\n");
    return
end

% Plot graphs
switch graph
    case 'sep'
        figure;
        %plot(t_control,y_control(1,:),'r--')
        %hold on
        plot(t,y(1,:),'b-')
        %plot(t2,y2(1,:))
        xlabel('Time - $\tau$','Interpreter','latex','fontsize',14);
        ylabel('Valve position - $y_1$','Interpreter','latex','fontsize',14)
        %legend('Initial','Quarter Wave Model','location','northeast')
        %pbaspect([3 1 1])
        %set(gcf,'Position',[216 93 560 420])
        %set(gcf,'Position',[216 93 670 260])

        figure;
        %plot(t_control,y_control(2,:),'r--')
        %hold on
        plot(t,y(2,:),'b-')
        %plot(t2,y2(2,:))
        xlabel('Time'); ylabel('Valve velocity')
        %legend('Initial','Quarter Wave Model','location','northeast')

        figure;
        %plot(t_control,y_control(4,:),'r--')
        %hold on
        plot(t,y(4,:),'b-')
        %plot(t2,y2(4,:))
        xlabel('Time'); ylabel('Pressure fluctuation')
        %legend('Initial','Quarter Wave Model','location','northwest')

        figure;
        %plot(t_control,y_control(5,:),'r--')
        %hold on
        plot(t,y(5,:),'b-')
        %plot(t2,y2(5,:))
        xlabel('Time'); ylabel('Velocity fluctuation')
        %legend('Initial','Quarter Wave Model','location','northwest')
    case 'all'
        figure;
        set(gcf,'Position',[216 93 560 680])

        fig1 = subplot(4,1,1);
        %plot(t_control,y_control(1,:),'r--')
        %hold on
        plot(t,y(1,:),'b-')
        %plot(t2,y2(1,:))
        y1 = ylabel('Valve position');
        fig1.XTickLabel = [];

        fig2 = subplot(4,1,2);
        %plot(t_control,y_control(2,:),'r--')
        %hold on
        plot(t,y(2,:),'b-')
        %plot(t2,y2(2,:))
        y2 = ylabel('Valve velocity');
        fig2.XTickLabel = [];

        fig3 = subplot(4,1,3);
        %plot(t_control,y_control(4,:),'r--')
        %hold on
        plot(t,y(4,:),'b-')
        %plot(t2,y2(4,:))
        y3 = ylabel('Pressure fluctuation');
        fig3.XTickLabel = [];

        subplot(4,1,4)
        %plot(t_control,y_control(5,:),'r--')
        %hold on
        plot(t,y(5,:),'b-')
        %plot(t2,y2(5,:))
        xlabel('Time'); y4 = ylabel('Velocity fluctuation');
        %legend('Initial','Quarter Wave Model','location','northwest')

        % 
        yc1 = y1.Position; yc2 = y2.Position; yc3 = y3.Position; yc4 = y4.Position;
        xpos = min([yc1(1), yc2(1), yc3(1), yc4(1)]);
        yc1(1) = xpos; set(y1,'Pos',yc1)
        yc2(1) = xpos; set(y2,'Pos',yc2)
        yc3(1) = xpos; set(y3,'Pos',yc3)
        yc4(1) = xpos; set(y4,'Pos',yc4)
end

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