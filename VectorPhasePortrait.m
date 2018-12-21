%%
%% Initialise variables
a = 0.01;
b = 10;

%r = 0.9;

x_max = 1.5;
x_min = 0;

%% 
grid_points = [41,11,41];

x_range = linspace(x_min,x_max,grid_points(1));
v_range = linspace(-0.3,0.3,grid_points(2));
p_range = linspace(0,2,grid_points(3));

%% 
[X,V,P] = meshgrid(x_range, v_range, p_range);
DDT = Derivative(0,[X(:)';V(:)';P(:)'],a,b);
DXDT = DDT(1,:); DVDT = DDT(2,:); DPDT = DDT(3,:); clear DDT;
%DXDT = reshape(DXDT,grid_points);
%DVDT = reshape(DVDT,grid_points);
%DPDT = reshape(DPDT,grid_points);

%%
figure
quiver3(X(:)',V(:)',P(:)',DXDT,DVDT,DPDT)
hold on
plot3(X(:)',zeros(size(X(:)')),(1./(X(:).^2))','r--')
xlim([x_range(1),x_range(end)]); xlabel('Position')
ylim([v_range(1),v_range(end)]); ylabel('Velocity')
zlim([p_range(1),p_range(end)]); zlabel('Pressure')

%%
[X,V,P] = meshgrid(x_range, 0, p_range);
DDT = Derivative(0,[X(:)';V(:)';P(:)'],a,b);
DXDT = DDT(1,:); DVDT = DDT(2,:); DPDT = DDT(3,:); clear DDT;
DXDT = reshape(DXDT,[grid_points(1),grid_points(3)]);
%DVDT = reshape(DVDT,grid_points); % All equal to zero
DPDT = reshape(DPDT,[grid_points(1),grid_points(3)]);

size(X)
size(V)
size(P)

%%
figure
quiver(X(:),P(:),DXDT(:),DPDT(:))
hold on
x = linspace(x_range(1),x_range(end),101);
plot(x,1./(x.^2),'r--')
xlim([x_range(1),x_range(end)]); ylim([0,2])
xlabel('Position'); ylabel('Pressure')