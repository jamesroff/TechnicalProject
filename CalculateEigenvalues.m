function [eigen] = CalculateEigenvalues(pres,gamma,alpha,delta,sigma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%     %[gamma, q, Lambda, alpha, delta, kappa, beta, mu, sigma, phi] = p;
%     gamma = parameters(1);
%     q = parameters(2);
%     % Lambda = 0;
%     alpha = parameters(4);
%     delta = parameters(5);
%     % kappa = 0;
%     % beta = 0;
%     %mu = parameters(8);
%     sigma = parameters(9); 
%     % phi = 0;
    
    % Set characteristic polynomial
    %p = [1, ((pi*sqrt(2)*alpha*sigma)/(4*gamma*sqrt(pres))), ( ((pi^2*sqrt(2))/(8*gamma^2)) - 2*sqrt(2)*delta^2*pres + 2*pres ), ((3*pi*sqrt(2)*alpha*sigma*sqrt(pres))/(2*gamma)), ((pi^2*pres)/(2*gamma^2))];
    %eig = roots(p);
    
    % Jacobian
    k = 0;
    J = [0, 1, 0, 0;
         2*pres, -k, 2, 0;
         0, sqrt(2)*(delta^2)*pres, 0, (pi/2)*(alpha/gamma);
         0, -sigma*sqrt(pres)*(sqrt(2)+delta^2), -(pi/2)*(1/(alpha*gamma)), -((pi*sqrt(2))/4)*((alpha*sigma)/(gamma*sqrt(pres)))]
    [~,D] = eig(J);
    eigen = zeros(1,4);
    for i = 1:4
        eigen(i) = D(i,i);
    end
    eigen
     
    % Plot eigenvalues
    figure; hold on; grid on;
    for i = 1:4
        plot(real(eigen(i)),imag(eigen(i)),'kx')
    end
end

