function [A B O D B_plus B_minus]=ssmod_Galpha(phi,omega,gamma,Sa,Sd,Sp,dt,varargin)

%% State space model from G-alpha discretization
%
% From Aucejo (2019)
%
% Inputs:
% phi: mode shape matrix
% omega: natural frequency matrix with diagonal omega
% gamma: damping matrix with diagonal 2*omega*xi
% Sa: selection matrix for acceleration
% Sd: selection matrix for displacement
% Sp: selection matrix for input
% dt: time discretization
%
% Outputs:
% A: state matrix
% B: input matrix
% O: output matrix
% D: direct transmission matrix
% B_plus: other system matrix, see below
% B_minus: other system matrix, see below
%
% Model:
% x(k+1)=A*x(k)+Bplus*u(k+1)+Bminus*u(k)+w(k)
% y(k)=O*x(k)+v(k)
%
% Modal state: x(k)=[z(k) z(k)' z(k)'']^T
%
% Transformation:
% x(k)=xbar(k)+Bplus*u(k)
%
% Introduced to model:
% xbar(k+1)=A*xbar(k)+B*u(k)+w(k)
% y(k)=O*xbar(k)+D*u(k)+v(k)
% where
% B=A*B_plus+B_minus
% D=O*B_plus
%
%%

p=inputParser;
addParameter(p,'force','modal',@ischar)
parse(p,varargin{:});
force = p.Results.force;

%%

nm=size(omega,1);
Sv=sparse(size(Sa,1),size(Sa,2));
Z=gamma;

%%

p_inf=1;
alpha_f=p_inf/(1+p_inf);
alpha_m=3*alpha_f-1;
gamma=0.5+alpha_f-alpha_m;
beta=0.25*(1+alpha_f-alpha_m)^2;
h=dt;

Ld=eye(nm)/(    (1-alpha_m)/(beta*h^2)*eye(nm)	+	(1-alpha_f)*gamma/(beta*h)*Z	+   (1-alpha_f)*omega^2);

Kd=Ld*(    (1-alpha_m)/(beta*h^2)*eye(nm)  +   (1-alpha_f)*gamma/(beta*h)*Z	+   -alpha_f*omega^2);

Cd=Ld*(    (1-alpha_m)/(beta*h)*eye(nm)  +   (1-alpha_f)*(gamma/beta-1)*Z	+   -alpha_f*Z);

Md=Ld*(    (1-alpha_m)*(1/(2*beta)-1)*eye(nm)  +   (1-alpha_f)*(gamma/(2*beta)-1)*h*Z	+   -alpha_m*eye(nm));

Lv=gamma/(beta*h)*Ld;

Kv=gamma/(beta*h)*(Kd-eye(nm));

Cv=gamma/(beta*h)*Cd+(1-gamma/beta)*eye(nm);

Mv=gamma/(beta*h)*Md+(1-gamma/(2*beta))*h*eye(nm);

La=1/(beta*h^2)*Ld;

Ka=1/(beta*h^2)*(Kd-eye(nm));

Ca=1/(beta*h^2)*(Cd-h*eye(nm));

Ma=1/(beta*h^2)*Md-(1/(2*beta)-1)*eye(nm);

A=[Kd Cd Md ; Kv Cv Mv ; Ka Ca Ma];

if strcmpi(force,'modal')
	TempTermForce=eye(nm);
else
	TempTermForce=phi.'*Sp;
end

B_plus=(1-alpha_f)*[Ld;Lv;La]*TempTermForce; %*phi.'*Sp;

B_minus=alpha_f*[Ld;Lv;La]*TempTermForce; %*phi.'*Sp;

O=[Sd*phi Sv*phi Sa*phi];

B=A*B_plus+B_minus;
D=O*B_plus;

