function [dAc_approx,dBc_approx,dGc_approx,dJc_approx,dAd_approx,dBd_approx,dGd_approx,dJd_approx]...
          =linearized_delta_modal(Omega0,Xi0,Phi0,dOmega,dXi,dPhi,Sd,Sv,Sa,Sp,dt)

%%

np=size(Sp,2);
ny=size(Sa,1);
ns=size(Omega0,1)*2;

%%

if isvector(Omega0); Omega0=diag(Omega0); end
if isvector(Xi0); Xi0=diag(Xi0); end
if isvector(dOmega); dOmega=diag(dOmega); end
if isvector(dXi); dXi=diag(dXi); end

% Mg0=eye(length(Omega0))

if isempty(dOmega)
    dOmega=zeros(size(Omega0));
end

if isempty(dXi)
    dXi=zeros(size(Xi0));
end

if isempty(dPhi)
    dPhi=zeros(size(Phi));
end

if isempty(Sv)
    Sv=zeros(size(Sa));
end

%%

% Continuous
dAc_approx=[zeros(ns/2) zeros(ns/2) ; -2*Omega0*dOmega -2*(Omega0*dXi+Xi0*dOmega) ];

dBc_approx=[zeros(ns/2,np) ; dPhi.'*Sp ];

dGc_approx=[ Sd*dPhi-Sa*(Phi0*2*Omega0*dOmega+dPhi*Omega0^2)    ,    Sv*dPhi-Sa*2*(Phi0*Omega0*dXi+Phi0*dOmega*Xi0+dPhi*Omega0*Xi0) ];
    
dJc_approx=Sa*(Phi0*dPhi.'+dPhi*Phi0.')*Sp;

% Discrete
dAd_approx=[ 0.5*dt^2*eye(ns/2) ; dt*eye(ns/2) ]*[-2*Omega0*dOmega -(2*Omega0*dXi+2*dOmega*Xi0)];

dBd_approx=[ 0.5*dt^2*eye(ns/2) ; dt*eye(ns/2) ]*[dPhi.'*Sp];

dGd_approx=dGc_approx;

dJd_approx=dJc_approx;


