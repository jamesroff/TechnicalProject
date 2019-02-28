% % syms d q mu sigma a pi root2;
% % syms omega;
% % 
% % syms ratio;
% % %ratio = (q*d)/(mu*sigma);
% % 
% % coef3 = ((root2 * a)/(2));
% % coef2 = 0;
% % coef1 = a*ratio^2;
% % coef0 = -ratio^2;
% % 
% % f = coef3*omega^3 + coef2*omega^2 + coef1*omega + coef0;
% % f = coef0*omega^3 + coef1*omega^2 + coef2*omega + coef3;
% % 
% % %F = factor(coef3*omega^3 + coef2*omega^2 + coef1*omega + coef0)
% % %F = factor(f,omega)
% % 
% % %R = root(f,omega);
% % R = solve(f,omega,'MaxDegree',3)

syms a c d;
b = 0;
syms w;

f = a*w^3 + b*w^2 + c*w + d;
%f = a + b*w + c*w^2 + d*w^3;

R = solve(f,w,'MaxDegree',3);

simplify(R(1))

%R(1)