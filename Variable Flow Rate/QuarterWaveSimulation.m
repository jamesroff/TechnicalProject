function [t,y,t_col,t_end] = QuarterWaveSimulation(p,init,x_range,t_range,model,cont)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% p = [gamma, q, Lambda, alpha, delta, kappa, beta, mu, sigma, phi, r]
r = p(end);

t = t_range(1); y = init;

opts = odeset('Events',@(t,y) Collision(t,y,x_range) );

t_col = [];

while t(end) ~= t_range(2)
    
    init = y(end,:);
    init(2) = - r * init(2);
    
    if (init(1) - x_range(1) <= 1e-5) && (init(2) <= 1e-3)
        break   
    elseif (abs (init(1) - x_range(2)) <= 1e-5) && (init(2) >= -1e-3)
        break
    end
    
    switch model
        case 'Initial'
            [t1,y1] = ode45(@(t,y) Derivative(t,y,p(1:end-1)), [t(end),t_range(2)], init, opts);
        case 'QWM'
            [t1,y1] = ode45(@(t,y) FullDerivative(t,y,p(1:end-1)), [t(end),t_range(2)], init, opts);
        case 'FullQWM'
            [t1,y1] = ode45(@(t,y) ComplexDerivative(t,y,p(1:end-1)), [t(end),t_range(2)], init, opts);
        case 'FixedPressure'
            [t1,y1] = ode45(@(t,y) AltDerivative(t,y,p(1:end-1),2), [t(end),t_range(2)], init, opts);
    end
    
    t_col = [t_col, t1(end)];
        
    %fprintf('Restarting\n')
    
    t = [t; NaN; t1]; y = [y; NaN([1,size(y,2)]); y1];
    
    if cont == false
        break
    end
end

if t_col(end) == t_range(2)
    t_col = t_col(1:end-1);
end

t_end = CalculateEndTime(t_col);

y = y';

end

%% Functions

function [v,t,d] = Collision(~,y,x_range)
    v = [real( y(1) ) - x_range(1),...
         real( y(1) ) - x_range(2)];
    t = [1,1];
    d = [-1,+1];
end

function [t_end] = CalculateEndTime(t)
    if size(t,2) < 3
        t_end = Inf;
        return
    end
    r = (t(end)-t(end-1))/(t(end-1)-t(end-2));
    a = t(end-1) - t(end-2);
    timeToEnd = a/(1-r);
    t_end = t(end-2) + timeToEnd;
end