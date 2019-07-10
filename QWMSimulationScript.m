% Reset
close all;
%clear all;
%clc;

model = 'QWM'; % 'Initial' or 'QWM'
modes = 1;
valve = 'pilot';

% Parameters to change
cont = false; % set to true to simulate after initial collision

q = 0.6;
%%% Very weird behaviour with L = 9 - L = 10 - L = 11 for 'spring'
L = 20; % 0.3, 0.35, 0.5
x_range = [0,1.1];
t_range = [0,700];
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
pres = ((sqrt(2)*pi/4)*(p(4)/p(1))*(pi/(2*p(1)))^2)/(p(5)^2 * 2 * (1-a)); %+ 0.07;

init = zeros(1,3+2*modes);
init(1) = 1;
init(3) = pres;
init(4) = -0.001;

% Calculate QWM
[t_control, y_control] = QuarterWaveSimulation(modes,p,init,x_range,t_range,'Initial',true);
[t,y] = QuarterWaveSimulation(modes,p,init,x_range,t_range,model,cont);
%[t2,y2] = QuarterWaveSimulation(p,init,x_range,t_range,model);
if isempty(t) && isempty(y)
    fprintf("Error during Quarter Wave Simulation\n");
    return
end

% Create legend labels
Labels = {'Valve position','Valve velocity','Tank Pressure'};
for i = 1:modes
    Labels{3+i} = ['Pressure - mode ' num2str(i,'%d')];
    Labels{3+modes+i} = ['Velocity - mode ' num2str(i,'%d')];
end

% Plot graphs
for i = 1:size(y_control,1)
    figure;
    plot(t_control,y_control(i,:),'r--')
    hold on
    plot(t,y(i,:),'b-')
    %plot(t2,y2(i,:))
    xlabel('Time'); ylabel(Labels(i))
    
    if i < 4
        legend('Initial','Quarter Wave Model','location','northeast')
    else
        legend('Initial','Quarter Wave Model','location','northwest')
    end
    %pbaspect([3 1 1])
    %set(gcf,'Position',[216 93 560 420])
    %set(gcf,'Position',[216 93 670 260])
end