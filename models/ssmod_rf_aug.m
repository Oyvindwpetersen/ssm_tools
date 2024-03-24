function [Ad Bd Gd Jd Ac Bc Gc Jc F]=ssmod_rf_aug(K,C,M,Ahat,Bhat,Chat,Sa,Sd,Sp,dt)

%% State space model of augmented system with velocity input
%
% [udot ]   [    0            I            0     ][u    ]   [     0      ]
% [uddot] = [-inv(M)*K   -inv(M)*C   -inv(M)*Chat][udot ] + [  inv(M)*Sp ] p
% [sdot ]   [    0          Bhat        Ahat     ][s    ]   [     0      ]
%
% Inputs:
% K: stiffness matrix
% C: damping matrix
% M: mass matrix
% Ahat: rational function state matrix
% Bhat: rational function input matrix
% Chat: rational function output matrix
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
nrf=size(Ahat,1);

Ac=zeros(2*ndof+nrf,2*ndof+nrf);
Bc=zeros(2*ndof+nrf,np);
Jc=zeros(ny,2*ndof+nrf);
Gc=zeros(ny,np);

MinvK=M\K;
MinvC=M\C;
MinvChat=M\Chat;
Minv=eye(ndof)/M;

Ac=[zeros(ndof) eye(ndof) zeros(ndof,nrf) ;
    -M\K -M\C -M\Chat ;
    zeros(nrf,ndof) Bhat Ahat];
	
Bc=[zeros(ndof,np) ; Minv*Sp ; zeros(nrf,np)];

if isempty(Sd) & isempty(Sa)
	Gc=[];
	Jc=[];
else
	Gc=[(Sd-Sa*MinvK) , -Sa*MinvC -Sa*MinvChat ];
	Jc=[Sa * Minv * Sp];
end

[Ad Bd Gd Jd F]=ssmod_c2d(Ac,Bc,Gc,Jc,dt);

