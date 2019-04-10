%% Set up workspace
%close all;
%clear all;
%clc;

%% Initialise variables
a = 0.01;
b = 1;
r0 = [1,0,1];

r = 0.9;

t_range = [0,180];

x_max = 1.5; %x_max = 1000;
x_min = 0;

%%
x_range = [0,x_max];
%x_range = [0.4,1];
p_range = [0,10];
%p_range = [0,5];
v_range = [-0.3,0.3];

n = 5;

init = [linspace(x_range(1),x_range(2),n)', zeros(n,1), ones(n,1)*p_range(1);...
        ones(n,1)*x_range(2), zeros(n,1), linspace(p_range(1),p_range(2),n)';...
        linspace(x_range(2),x_range(1),n)', zeros(n,1), ones(n,1)*p_range(2)];
    
%init = [0.9,0,1];

%% Calculate trajectory
x = linspace(x_range(1),x_range(2),101);

figure;
plot(x,(1./(x.^2)),'r--')
hold on
for i = 1:size(init,1)
    [t,y] = InitialModel(a,b,r,init(i,:),[x_min,x_max],t_range);
    plot(y(1,:),y(3,:),'g-')
end
xlabel('Position - $\tilde{x}$','Interpreter','LaTeX');
ylabel('Pressure - $\tilde{p}$','Interpreter','LaTeX');
ylim(p_range); xlim(x_range)

%% 3D Init
x_values = linspace(x_range(1),x_range(2),2*n+1);
v_values = linspace(v_range(1),v_range(2),n);
p_values = linspace(p_range(1),p_range(2),n);

init = CreateInitialConditions(x_values,v_values,p_values);


%% Calculate trajectory
x = linspace(x_range(1),x_range(2),101);

figure;
plot3(x,zeros(size(x)),(1./(x.^2)),'r--')
hold on
for i = 1:size(init,1)
    [t,y] = InitialModel(a,b,r,init(i,:),[x_min,x_max],t_range);
    plot3(y(1,1),y(2,1),y(3,1),'go')
    plot3(y(1,:),y(2,:),y(3,:),'g-')
end
xlabel('Position - $\tilde{x}$','Interpreter','LaTeX');
ylabel('Velocity - $\tilde{x}''$','Interpreter','LaTeX');
zlabel('Pressure - $\tilde{p}$','Interpreter','LaTeX');
xlim(x_range); ylim(v_range); zlim(p_range)

function [init] = CreateInitialConditions(x_vals,v_vals,p_vals)
    %init = [0.99, 0, 1];
    x = length(x_vals); v = length(v_vals); p = length(p_vals);
    init = zeros(v*p+x*v+2*x*p,3);
    
    % Bottom surface - p = 0
    [X,V] = meshgrid(x_vals,v_vals);
    init(1:x*v,:) = [X(:), V(:), p_vals(1)*ones(size(X(:)))];
    
    % Side surface 1 - v = -0.3
    [X,P] = meshgrid(x_vals,p_vals);
    init(x*v+1:x*v+x*p,:) = [X(:), v_vals(1)*ones(size(X(:))), P(:)];
    
    % Side surface 2 - v = 0.3
    init(x*v+x*p+1:x*v+2*x*p,:) = [X(:), v_vals(end)*ones(size(X(:))), P(:)];
    
    % Back surface - x = 1
    [V,P] = meshgrid(v_vals,p_vals);
    init(x*v+2*x*p+1:v*p+x*v+2*x*p,:) = [x_vals(end)*ones(size(V(:))), V(:), P(:)];
end