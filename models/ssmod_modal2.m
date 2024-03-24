function [Ad Bd Gd Jd Ac Bc Gc Jc F]=ssmod_modal2(Mg,phi,omega,gamma,Sa,Sd,Sp,dt,varargin)

%% Impulse response function for state space model

% Inputs:
% Mg: generalized mass
% phi: mode shape matrix
% omega: natural frequency matrix with diagonal omega
% gamma: damping matrix with diagonal 2*omega*xi
% Sa: selection matrix for acceleration
% Sd: selection matrix for displacement
% Sp: selection matrix for input
% dt: time discretization

% Outputs:
% Ad: state matrix (disc)
% Bd: input matrix (disc)
% Gd: output matrix (disc)
% Jd: direct transmission matrix (disc)
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
	TempTermForce=Mg\eye(nm);
    np=nm;
elseif strcmpi(force,'disc')
	TempTermForce=Mg\phi.'*Sp;
else
    error('Missing arguement (modal or disc)')
end

Bc=[zeros(nm,np);TempTermForce];

Gc=[Sd*phi-Sa*phi*omega.^2, Sv*phi-Sa*phi*gamma];
Jc=[Sa*phi*TempTermForce];

Ad=expm(dt*Ac);
Bd=[Ad-eye(size(Ad))]*(Ac\Bc);

Gd=Gc;
Jd=Jc;

%first order hold
F=Ac\(Bd-Bc*dt)*dt^-1;
