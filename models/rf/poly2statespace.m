function [Ac,Bc,Cc,Dc]=poly2statespace(p,q)

%% State space model from poly rational function
%
% Inputs:
% p: polynomial coefficients
% q: polynomial coefficients
%
% q=[q_n q_(n-1) ... q_1 q_0] % n+1 terms 
% p=[p_r p_(r-1) ... p_1 p_0] % r+1 terms 
%
% Outputs:
% Ac: state matrix (cont)
% Bc: input matrix (cont)
% Cc: output matrix (cont)
% Dc: direct transmission matrix (cont)
%
%%

n=length(q)-1; % highest poly order in q
r=length(p)-1; % highest poly order in p

Ac=[ zeros(n-1,1) eye(n-1) ; -flip(q(2:end)) ];
Bc=zeros(n,1); Bc(end)=1./q(1);

Cc=zeros(1,n);
Cc(1,1:(r+1))=flip(p);
Dc=0;