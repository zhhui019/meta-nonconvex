function y = etaprimemcp(x, lambda, a)
% etaprimemcp is the derivative of the mcp threshold function.

y = (abs(x) > a*lambda) + (a/(a-1))*(abs(x) <= a*lambda & abs(x) > lambda) ;
