function [] = CalculateFrequency(gamma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

w = (pi)/(2*gamma); % w = 2 pi f
fprintf('\nThe quarter-wave angular frequency is w = %d\n',w)
f = w/(2*pi);
fprintf('The quarter-wave frequency is f = %d\n',f)
fprintf('The quarter-wave period is %d\n\n',1/f)

end