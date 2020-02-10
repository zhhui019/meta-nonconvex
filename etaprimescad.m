function y = etaprimescad(x, lambda)
% etaprimescad is the derivative of the scad threshold function. 
a=3.7;

y = (abs(x) > a*lambda) + (a-1)/(a-2)*(abs(x) <= a*lambda & abs(x) > 2*lambda) + (abs(x) <= 2*lambda & abs(x) > lambda);