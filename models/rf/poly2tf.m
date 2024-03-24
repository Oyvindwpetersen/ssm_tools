function H=poly2tf(omega_axis,p,q)

%% Transfer function for poly rational function
%
% Inputs:
% omega_axis: frequency axis in [rad/s]
% p: polynomial coefficients
% q: polynomial coefficients
%
% Outputs:
% H: transfer function
%
%%

P=polyval(p,1i*omega_axis);

Q=polyval(q,1i*omega_axis);

H(1,1,:)=P./Q;