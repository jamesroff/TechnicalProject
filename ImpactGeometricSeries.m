% 
close all;
clear all;
clc;

% Calculate valve motion
p = [0.2212, 0.6, 0, 8.566, 1, 0, 0.0433, 0.141, 10.38, 0, 0.8];
pres = ((p(2))/(p(8)*(9)))^2;
[t,y,t_col,t_end3] = QuarterWaveSimulation(p,[1.01,0,pres,0,0],[0,1.5],[0,25],'Initial');

% Geometric series calculations
geo = t_col(2:end)-t_col(1:end-1);
ratio = geo(2:end)./geo(1:end-1);

r = (t_col(3)-t_col(2))/(t_col(2)-t_col(1));
a = t_col(2) - t_col(1);
t_end1 = t_col(1) + a/(1-r);

r = mean(ratio);
a = geo(1);
t_end2 = t_col(1) + a/(1-r);

% Valve position graph
figure;
plot(t,y(1,:),'k-')
hold on
y_lim = ylim;
plot([1,1]*t_end1,y_lim,'r--')
plot([1,1]*t_end2,y_lim,'b--')
plot([1,1]*t_end3,y_lim,'g--')
ylim(y_lim);

% geometric series
figure;
plot(ratio,'k-') % plot common ratio
hold on
x_lim = xlim;
plot(xlim,[1,1]*ratio(1),'r--')
plot(xlim,[1,1]*mean(ratio),'b--')
plot(xlim,[1,1]*ratio(end),'g--')
xlabel('Index of impact'); ylabel('Common ratio')

% time
figure;
plot(t_col,1:length(t_col),'k-')
hold on
plot([1,1]*t_end1,[1,length(t_col)],'r--')
plot([1,1]*t_end2,[1,length(t_col)],'b--')
plot([1,1]*t_end3,[1,length(t_col)],'g--')
xlabel('Time'); ylabel('Index of collision')