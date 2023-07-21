function [Fc,Lc,Hc,sigma_w]=ssmod_matern(lambda,p_order,sigma_p)

%% State space model for Matern
%
% ds/dt = F*s(t) + L*w(t)
% p(t) = H*s(t)
%
% Inputs:
% lambda: hyperparameter, inverse length scale
% p_order: order of the differential equation (+1), v=p+0.5
% sigma_p: standard deviation of p
%
% Outputs:
% Fc: state matrix (cont)
% Lc: input matrix (cont)
% Hc: output matrix (cont)
%
% Depending on the order of p:
% p=0: CF(tau) ~ c0*exp(-b*tau) (Ornstein-Uhlenbeck)
% p=1: CF(tau) ~ (c0+c1*tau)*exp(-b*tau)
% p=2: CF(tau) ~ (c0+c1*tau+c2*tau^2)*exp(-b*tau)
%

%% Assign state space matrices

if p_order==0
	Fc=-lambda;
elseif p_order==1
	Fc=[0 1 ; -lambda.^2 -2*lambda];
elseif p_order==2
	Fc=[0 1 0 ; 0 0 1 ; -lambda.^3 -3*lambda^2 -3*lambda];
elseif p_order==3
	Fc=[0 1 0 0 ; 0 0 1 0 ; 0 0 0 1 ; -lambda.^4 -4*lambda^3 -6*lambda^2 -4*lambda];
else
	error('p must any of 0,1,2,3')
end

Lc=zeros(p_order+1,1); Lc(end)=1;
Hc=zeros(1,p_order+1); Hc(1)=1;


%% Variance of p:

% From Maple:

% sigma_p^2=sigma_w^2 * 1/(2*lambda)
% sigma_p^2=sigma_w^2 * 1/(4*lambda^3)
% sigma_p^2=sigma_w^2 * 3/(16*lambda^5)
% sigma_p^2=sigma_w^2 * 5/(32*lambda^7)

if p_order==0
	sigma_w=sqrt( sigma_p.^2 * (2*lambda) );
elseif p_order==1
	sigma_w=sqrt( sigma_p.^2 * (4*lambda^3) );
elseif p_order==2
	sigma_w=sqrt( sigma_p.^2 * 1/3*(16*lambda^5) );
elseif p_order==3
	sigma_w=sqrt( sigma_p.^2 * 1/5*(32*lambda^7) );
else
	error('p must any of 0,1,2,3')
end