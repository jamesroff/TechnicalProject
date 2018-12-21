function [dydt] = Derivative(~,y,a,b)
    dydt = zeros(size(y));
    
    dydt(1,:) = y(2,:);
    dydt(2,:) = - y(2,:) + a.*y(3,:).*( y(1,:).^2 - 1 );
    dydt(3,:) = b*(1 - y(1,:).*sqrt(y(3,:)));
end