function [dydt] = Derivative(~,y,p)
    dydt = zeros(size(y));
    
    %[gamma, q, Lambda, alpha, delta, kappa, beta, mu, sigma, phi] = p;
    %gamma = p(1);
    q = p(2);
    %Lambda = p(3);
    %alpha = p(4);
    %delta = p(5);
    kappa = p(6);
    beta = p(7);
    mu = p(8);
    sigma = p(9);
    %phi = p(10);
    
    y(4,:) = 0;
    y(5,:) = 0;
    
    dydt(1,:) = y(2,:);
    dydt(2,:) = - kappa * y(2,:) + y(3,:).*(y(1,:).^2 - 1) + 2*y(4,:).*y(1,:).^2;
    dydt(3,:) = beta*( q - sigma*y(5,:) - sigma*mu*y(1,:).*sqrt(y(3,:) + y(4,:)) );
    
    dydt(4,:) = 0;
    dydt(5,:) = 0;
end