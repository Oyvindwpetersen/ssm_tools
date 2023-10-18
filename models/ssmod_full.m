function [A B G J Ac Bc Gc Jc F]=ssmod_full(K,C,M,Sa,Sd,Sp,dt)

%% State space model of full (non-reduced) (M,C,K) system
%
% Inputs:
% K: stiffness matrix
% C: damping matrix
% M: mass matrix
% Sa: selection matrix for acceleration
% Sd: selection matrix for displacement
% Sp: selection matrix for input
% dt: time discretization
%
% Outputs:
% A: state matrix (disc)
% B: input matrix (disc)
% G: output matrix (disc)
% J: direct transmission matrix (disc)
% Ac: state matrix (cont)
% Bc: input matrix (cont)
% F: matrix for first-order-hold
%
%% 

ndof=size(K,1);
ny=size(Sa,1);
np=size(Sp,2);

Ac=zeros(2*ndof,2*ndof);
Bc=zeros(2*ndof,np);
Jc=zeros(ny,2*ndof);
Gc=zeros(ny,np);

MinvK=M\K;
MinvC=M\C;
Minv=eye(ndof)/M;

Ac=[zeros(ndof) eye(ndof) ; -MinvK -MinvC];

if isempty(Sd) & isempty(Sa)
    
else
    Bc=[zeros(ndof,np) ; Minv*Sp];
end


if isempty(Sd) & isempty(Sa)
	Gc=[];
	Jc=[];
else
	Gc=[(Sd-Sa*MinvK) , -Sa*MinvC ];
	Jc=[Sa * Minv * Sp];
end

[A B G J F]=ssmod_c2d(Ac,Bc,Gc,Jc,dt);

