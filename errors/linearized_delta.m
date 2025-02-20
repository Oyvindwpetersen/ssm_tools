function [dAc_approx,dBc_approx,dGc_approx,dJc_approx,dAd_approx,dBd_approx,dGd_approx,dJd_approx]=linearized_delta(M0,C0,K0,dM,dC,dK,Sa,Sp,dt)
%% Approximation deltas in state-space matrices due to perturbations in mass, damping, or stiffness
%
% Inputs:
% M0: erroneous mass matrix
% C0: erroneous damping matrix
% K0: erroneous stiffness matrix
% dM: change in mass matrix
% dC: change in damping matrix
% dK: change in mstiffnessass matrix
% Sa: selection matrix for accelerations
% Sp: selection matrix for forces
% dt: time step
% disc_type: 'zoh' or 'foh'
%
% Outputs:
% dAc_approx: change in state matrix (cont)
% dBc_approx: change in input matrix (cont)
% dGc_approx: change in output matrix (cont)
% dJc_approx: change in direct transmission matrix (cont)
% dAd_approx: change in state matrix (disc)
% dBd_approx: change in input matrix (disc)
% dGd_approx: change in output matrix (disc)
% dJd_approx: change in direct transmission matrix (disc)
%
%%

n_dof=size(M0,1);
ny=size(Sa,1);
ns=size(M0,1)*2;
np=size(Sp,2);

%%

if isempty(dM)
    dM=zeros(size(M0));
end

if isempty(dC)
    dC=zeros(size(C0));
end

if isempty(dK)
    dK=zeros(size(K0));
end

%%

% Continuous
dAc_approx=[zeros(n_dof) zeros(n_dof) ; (-M0\dK+M0\dM/M0*K0) (-M0\dC+M0\dM/M0*C0) ];

dBc_approx=[zeros(n_dof,np) ; (-M0\dM/M0)*Sp ];

dGc_approx=[Sa*(-M0\dK+(M0\dM/M0)*K0)   Sa*(-M0\dC+(M0\dM/M0)*C0) ];

dJc_approx=Sa*(-M0\dM/M0)*Sp;

% Discrete
dAd_approx=[ 0.5*dt^2*eye(ns/2) ; dt*eye(ns/2) ]*[(-M0\dK+M0\dM/M0*K0) (-M0\dC+M0\dM/M0*C0) ];

dBd_approx=[ 0.5*dt^2*eye(ns/2) ; dt*eye(ns/2) ]*[ (-M0\dM/M0)*Sp];

dGd_approx=dGc_approx;

dJd_approx=dJc_approx;



% dAd_approx=dt*[dt/2*(-M0\dK+M0\dM/M0*K0) dt/2*(-M0\dC+M0\dM/M0*C0) ; (-M0\dK+M0\dM/M0*K0) (-M0\dC+M0\dM/M0*C0) ];
% 
% dBd_approx=dt*[-dt/2*(M0\dM/M0)*Sp ; -(M0\dM/M0)*Sp ];
% 
% dGd_approx=dGc_approx;
% 
% dJd_approx=dJc_approx;


% if strcmpi(disc_type,'zoh')
%     dBd_approx=dt*[-dt/2*(M0\dM/M0)*Sp ; -(M0\dM/M0)*Sp ];
%     dJd_approx=dJc_approx;
% elseif strcmpi(disc_type,'foh')
%     dBd_approx=dt*[-dt*(M0\dM/M0)*Sp ; -(M0\dM/M0)*Sp ];
%     dJd_approx=dJc_approx;
% end

    



