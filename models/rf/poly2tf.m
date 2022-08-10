function H=poly2tf(w_axis,p,q)

%% Transfer function for poly rational function
%
% Inputs:
% w_axis: frequency axis in [rad/s]
% p: polynomial coefficients
% q: polynomial coefficients
%
% Outputs:
% H: transfer function
%
%%

P=polyval(p,1i*w_axis);

Q=polyval(q,1i*w_axis);

H(1,1,:)=P./Q;