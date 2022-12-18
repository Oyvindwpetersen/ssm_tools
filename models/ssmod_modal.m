function [A B G J Ac Bc Gc Jc F]=ssmod_modal(phi,omega,gamma,Sa,Sd,Sp,dt,varargin)

%% Impulse response function for state space model

% Inputs:
% phi: mode shape matrix
% omega: natural frequency matrix with diagonal omega
% gamma: damping matrix with diagonal 2*omega*xi
% Sa: selection matrix for acceleration
% Sd: selection matrix for displacement
% Sp: selection matrix for input
% dt: time discretization

% Outputs:
% A: state matrix (disc)
% B: input matrix (disc)
% G: output matrix (disc)
% J: direct transmission matrix (disc)
% Ac: state matrix (cont)
% Bc: input matrix (cont)
% F: matrix for first-order-hold

%%

p=inputParser;
addParameter(p,'force','modal',@ischar)
parse(p,varargin{:});
force = p.Results.force;

%%

nm=size(omega,1);
ny=size(Sa,1);
np=size(Sp,2);

Ac=zeros(2*nm,2*nm);
Bc=zeros(2*nm,np);
Jc=zeros(ny,2*nm);
Gc=zeros(ny,np);

Sv=sparse(ny,size(phi,1));

% Reduced order state space model
Ac=[zeros(nm,nm) eye(nm); -omega.^2 -gamma];

if strcmpi(force,'modal')
	TempTermForce=eye(nm);
    np=nm;
elseif strcmpi(force,'disc')
	TempTermForce=phi.'*Sp;
else
    error('Missing arguement (modal or disc)')
end

Bc=[zeros(nm,np);TempTermForce];

Gc=[Sd*phi-Sa*phi*omega.^2, Sv*phi-Sa*phi*gamma];
Jc=[Sa*phi*TempTermForce];

A=expm(dt*Ac);
B=[A-eye(size(A))]*(Ac\Bc);

G=Gc;
J=Jc;

%first order hold
F=Ac\(B-Bc*dt)*dt^-1;
