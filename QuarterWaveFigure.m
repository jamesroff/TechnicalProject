%% Initialise coordinates
L = 1;
x = linspace(0,L,1001);

p_0 = 1;
v_L = 1;

B_1 = 0.2;
B_3 = 0.1;

C_1 = 0.2;
C_3 = 0.1;

%% Plot figure
figure;
fig1 = subplot(2,1,1);
plot(x,p_0*ones(size(x)),'k-')
hold on; grid on
plot(x,p_0 + B_1*sin((pi/(2*L))*x),'b--')
plot(x,p_0 + B_3*sin((pi/(2*L))*3*x),'r-.')
fig1.XTickLabel = []; fig1.YTickLabel = [];
ylim([0.85,1.25])

fig2 = subplot(2,1,2);
plot(x,v_L*ones(size(x)),'k-')
hold on; grid on
plot(x,p_0 + C_1*cos((pi/(2*L))*x),'b--')
plot(x,p_0 + C_3*cos((pi/(2*L))*3*x),'r-.')
fig2.XTickLabel = []; fig2.YTickLabel = [];
ylim([0.85,1.25])