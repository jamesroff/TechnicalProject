function [t,y] = QuarterWaveSimulation(p,init,x_range,t_range)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% p = [gamma, q, Lambda, alpha, delta, kappa, beta, mu, sigma, phi, r]
r = p(end);

t = t_range(1); y = init;

opts = odeset('Events',@(t,y) Collision(t,y,x_range) );

while t(end) ~= t_range(2)
    
    init = y(end,:);
    init(2) = - r * init(2);
    
    if (init(1) - x_range(1) <= 1e-5) && (init(2) <= 1e-3)
        break
    end

%     if init(1) - x_range < 0
%         break
%     end
    
    [t1,y1] = ode45(@(t,y) Derivative(t,y,p(1:end-1)), [t(end),t_range(2)], init, opts);
    
    %fprintf('Restarting\n')
    
    t = [t; NaN; t1]; y = [y; NaN([1,size(y,2)]); y1];
end

y = y';

end

%% Functions

% Derivative function
% function [dydt] = Derivative(~,y,a,b)
%     dydt = zeros(size(y));
%     size(y)
%     
%     dydt(1,:) = y(2,:);
%     dydt(2,:) = - y(2,:) + a.*y(3,:).*( y(1,:).^2 - 1 );
%     dydt(3,:) = b*(1 - y(1,:).*sqrt(y(3,:)));
% end

% Event function
% function [v,t,d] = Collision(~,y,x_range)
%     v = [y(1) - x_range(1),...
%          y(1) - x_range(2),...
%          x_range(1) < y(1) && y(1) < x_range(2)];
%     %y(1)
%     %x_range(1) < y(1) && y(1) < x_range(2)
%     t = [1,1,1];
%     d = [-1,+1,0];
% end

function [v,t,d] = Collision(~,y,x_range)
    v = [y(1) - x_range(1),...
         y(1) - x_range(2)];
    t = [1,1];
    d = [-1,+1];
end