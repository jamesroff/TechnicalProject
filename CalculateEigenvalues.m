function [eig] = CalculateEigenvalues(pres, parameters)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %[gamma, q, Lambda, alpha, delta, kappa, beta, mu, sigma, phi] = p;
    gamma = parameters(1);
    q = parameters(2);
    % Lambda = 0;
    alpha = parameters(4);
    delta = parameters(5);
    % kappa = 0;
    % beta = 0;
    %mu = parameters(8);
    sigma = parameters(9); 
    % phi = 0;
    
    % Set characteristic polynomial
    p = [1, ((pi*sqrt(2)*alpha*sigma)/(4*gamma*sqrt(p))), ( ((pi^2*sqrt(2))/(8*gamma^2)) - 2*sqrt(2)*delta^2*p + 2*p ), ((3*pi*sqrt(2)*alpha*sigma*sqrt(p))/(2*gamma)), ((pi^2*p)/(2*gamma^2))];
    eig = roots(p);

end

