function [dydt] = FullDerivative(~,y,p)
    dydt = zeros(size(y));
    
    %[gamma, q, Lambda, alpha, delta, kappa, beta, mu, sigma, phi] = p;
    gamma = p(1);
    q = p(2);
    Lambda = p(3);
    alpha = p(4);
    delta = p(5);
    kappa = p(6);
    beta = p(7);
    mu = p(8);
    sigma = p(9);
    phi = p(10);
    
    %pres = ((p(2))/(p(8)*p(9)))^2
    %y(3) = ((q)/(mu*sigma))^2;
    
    
    %root = sqrt( abs( y(3,:) + y(4,:) ) );
    root = sqrt( y(3,:) + y(4,:) );
    
    dydt(1,:) = y(2,:);
    dydt(2,:) = - kappa * y(2,:) + y(3,:).*(y(1,:).^2 - 1) + 2*y(4,:).*(y(1,:).^2);    
    
    dydt(3,:) = beta*( q - sigma*y(5,:) - sigma*mu*y(1,:).*root );
    
%     dydt(4,:) = (1/(1 - (sqrt(2)/2)*(delta^2)*(y(1,:).^2))) * ...
%                 ( (pi/2)*(alpha/gamma)*y(5,:)...
%                   - sqrt(2)*dydt(3,:)...
%                   + (sqrt(2)/2)*(delta^2)*(y(1,:).^2).*dydt(3,:)...
%                   + sqrt(2)*(delta^2)*y(1,:).*y(2,:).*( y(3,:) + y(4,:) )...
%                   + (pi/2)*((delta*Lambda)/(alpha))*y(1,:).*y(4,:).*root...
%                   - ((pi*sqrt(2))/(4))*Lambda*y(4,:).*y(5,:)                );
    dydt(4,:) = (pi/2)*(alpha/gamma)*y(5,:)...
                - sqrt(2)*dydt(3,:)...
                + (sqrt(2)/2)*(delta^2)*(y(1,:).^2).*dydt(3,:)...
                + sqrt(2)*(delta^2)*y(1,:).*y(2,:).*y(3,:)...
                + (pi/2)*((delta*Lambda)/(alpha))*y(1,:).*y(4,:).*root...
                - ((pi*sqrt(2))/(4))*Lambda*y(4,:).*y(5,:);
    
%     dydt(5,:) = - (pi/2)*(1/(alpha*gamma))*y(4,:)... %UNCHANGED
%                 + (pi/2)*Lambda*sigma * y(1,:).*y(5,:).*root ... %UNCHANGED
%                 + ((pi*sqrt(2))/(4))*Lambda* (y(5,:).^2) ... %UNCHANGED
%                 + (1/(1 - (sqrt(2)/2)*(delta^2)*(y(1,:).^2))) * ...
%                     ( - ((pi*sqrt(2))/(4)) * ((alpha*sigma)/(gamma)) * ((y(1,:).*y(5,:))./root)   ... % multiplied by B_dot term
%                     + (pi/4)*Lambda*sigma * (( y(1,:).*y(4,:).*y(5,:) )./root)                    ... % multiplied by B_dot term
%                     + ((pi*sqrt(2))/(4))*Lambda*(sigma^2) * (y(1,:).^2).*y(4,:)                   ... % multiplied by B_dot term
%                     - sqrt(2)*sigma* y(2,:).*root ...                                                 % NEW TERM
%                     + ((2 - sqrt(2))/(2))*sigma * ((y(1,:).*dydt(3,:))./(root))                 ) ... % NEW TERM
%                 + (sqrt(2)/2)*phi* (sqrt(2)*sigma*y(1,:).*root + y(5,:)) .* abs((sqrt(2)*sigma*y(1,:).*root + y(5,:))) ; %UNCHANGED
    dydt(5,:) = - (pi/2)*(1/(alpha*gamma))*y(4,:)...
                ... - ((pi*sqrt(2))/(4)) * ((alpha*sigma)/(gamma)) * ((y(1,:).*y(5,:))./root) ; ... THIS TERM IS THE PROBLEM!!!!!
                - sqrt(2)*sigma * y(2,:).*root ...
                + ((2 - sqrt(2))/2)*sigma*((y(1,:).*dydt(3,:))./root) ...
                - sigma*(delta^2) * (( (y(1,:).^2).*y(2,:).*y(3,:) )./root) ... % end of first line of eqn. 17
                - (1/2)*sigma*(delta^2) * (( (y(1,:).^3).*dydt(3,:) )./root) ...
                + (pi/2)*Lambda*sigma * y(1,:).*y(5,:).*root ...
                + ((pi*sqrt(2))/(4))*Lambda* (y(5,:).^2) ...
                + (pi/4)*Lambda*sigma * (( y(1,:).*y(4,:).*y(5,:) )./root) ...
                + ((pi*sqrt(2))/(4))*Lambda*(sigma^2) * (y(1,:).^2).*y(4,:) ... % end of second line of eqn. 17
                + (sqrt(2)/2)*phi* (sqrt(2)*sigma*y(1,:).*root + y(5,:)) .* abs((sqrt(2)*sigma*y(1,:).*root + y(5,:))) ... ; %THIS IS NON-ZERO
                ;
end