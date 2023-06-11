function y=polyvaleven(p,x)
%% Evaluate polynomial where all odd coefficients are zero
%
% y=p0+p2*x^2+p4*x^4...
%
% Inputs:
% p: polynomial coefficients for even terms, [p(k) p(k-2) ... p(2) p(0)]
% x: vector
%
% Outputs:
% y: polynomial
%
%% Evaluate polynomial recursively

% Copied from polyval.m, x replaced with x_squared 

x_squared=x.^2;

y=ones(size(x))*p(1);

nc=length(p);
for i=2:nc
    y=x_squared.*y + p(i);
end