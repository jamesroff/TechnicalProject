%% Set up workspace
%close all;
%clear all;
%clc;

%% Initialise variables
a = 0.01;
b = 1;
r0 = [1,0,1];
%init = r0 - [0.1,0.1,-0.1];
%init = [0,0,0];

r = 0.9;

t_range = [0,1000];

x_max = 1.5; x_max = 1000;
x_min = 0;

%% normal eigenvector for a = 0 and orthogonal directions
n = ((2*sqrt(5))/5)*[-0.5,0,1];
d1 = (sqrt(5)/5)*[2,0,1]; d2 = [0,1,0];
v1 = [0;0;1]; v2 = [(2-b)/(2*b);(b-2)/(2*b);1];

init = r0 + 0.1*v2';

%% Calculate trajectory
[t,y] = InitialModel(a,b,r,init,[x_min,x_max],t_range);

%% Plot graphs
figure;
plot(t,y(1,:))

figure;
plot(t,y(2,:))

figure;
plot(t,y(3,:))

%%
% figure
% plot3(y(1,:),y(2,:),y(3,:))
% hold on
% scatter(1,0,1,'rx')
% xlabel('Piston Position'); ylabel('Piston Velocity'); zlabel('Pressure')
% xlim([x_min,x_max]); %zlim([0,1])

%% 
orth = zeros(2,size(y,2));
   a = zeros(1,size(y,2));
for i = 1:length(y)
    %orth(:,i) = ProjectOrth(y(:,i),r0',n',d1',d2');
    [orth(:,i),a(i)] = ProjectEigen(y(:,i),r0',n',v1,v2);
end

figure
plot(orth(1,:),orth(2,:),'bx-')
%xlabel('Orthogonal Direction'); ylabel('Velocity')

figure
plot(t,a,'r-')

function [v] = ProjectOrth(R,r0,n,d1,d2)
    % vector given by :
    % R = r0 + a*n + x*d1 + y*d2
    a = dot(n,R-r0)
    x = dot(d1,R-r0); y = dot(d2,R-r0);
    v = [x;y];
end

function [v,a] = ProjectEigen(R,r0,n,v1,v2)
    % A d = B
    A = [ dot(n,n),  dot(n,v1),  dot(n,v2);
         dot(v1,n), dot(v1,v1), dot(v1,v2);
         dot(v2,n), dot(v2,v1), dot(v2,v2)];
    B = [dot(n,R-r0); dot(v1,R-r0); dot(v2,R-r0)];
    d = A\B;
    v = d(2:3);
    a = d(1);
end