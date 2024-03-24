function [Fc,Lc,Hc,sigma_w,S_trunc]=ssmod_squaredexp(L,sigma_p,ns,omega_axis)

%% State-space model approximation of SE kernel
%
% State equation:
% ds/dt=Fc*s(t)+Lc*w(t)
% p(t)=Hc*s(t)
%
% Inputs:
% L: length scale
% sigma_p: magnitude factor
% ns: state-space model order, N=ns+1 terms in Taylor series
% omega_axis: axis for spectral density
%
% Outputs:
% Fc: state kernel matrix in state-space model
% Lc: input matrix in state-space model
% Hc: output matrix in state-space model
% sigma_w: standard deviation of white noise input w(t) in state-space model 
% S_trunc: two-sided PSD of approximated model
%
% CF for SE model:
% R(t)=sigma_p^2*exp(-0.5*t^2/L^2)
%
% PSD from the inverse Fourier transform:
% S(w)=1/(2*pi) int R(t) exp(-i*w*tau) dtau
% S(w)=sigma_p^2*L/sqrt(2*pi)*exp(-0.5*L^2*w^2)
%
%% Checks

N=ns+1;

if ns>7 | ns<2
    error('ns must be in [2,7]');
end

if L<=0
	error('L must be greater than zero');
end
	
if sigma_p<=0
	error('sigma_p must be greater than zero');
end

%%

% Taylor expansion of exp(-0.5*L^2*w^2): 1/(w^0*L^0*A0+w^2*L^2*A2+w^4*L^4*A4+...until...w^10*L^10*A10)

% Powers and constants A for the polynomial
Power_even=[0 2 4 6 8 10 12 14];
A_coeff_even=[1 1/2 1/8 1/48 1/384 1/3840 1/46080 1/645120];
SignChange=[1 -1   1 -1   1 -1   1 -1];

% Truncate at N terms
Power_even_trunc=Power_even(1:N);
A_coeff_even_trunc=A_coeff_even(1:N);

% Rescale coeff so that highest poly has coeff unity
Abar_coeff_even_trunc=A_coeff_even_trunc.*L.^Power_even_trunc/(A_coeff_even_trunc(end)*L^Power_even_trunc(end));

% Constant value on top of fraction
Constant0=sigma_p^2*L/sqrt(2*pi)/(A_coeff_even_trunc(end)*L^Power_even_trunc(end));

% Fill zero coeff for odd w^1,w^3,...
Abar_coeff_all_trunc=[Abar_coeff_even_trunc .* SignChange(1:N) ; zeros(size(Abar_coeff_even_trunc))];
Abar_coeff_all_trunc=Abar_coeff_all_trunc(1:end-1);

% Roots of poly
roots_iw=roots(flip(Abar_coeff_all_trunc));

roots_iw_sort=cplxpair(roots_iw);
roots_iw_selected=roots_iw_sort(real(roots_iw_sort)<0);

poly_coeff=poly(roots_iw_selected);

%%

% The PSD of the white noise input is S_w=sigma_w^2/(2*pi)
% This must equal sigma_p^2*L/sqrt(2*pi)/(A_coeff_even_trunc(end)*L^Power_even_trunc(end))

sigma_w_squared=sqrt(2*pi)*sigma_p^2*L/(A_coeff_even_trunc(end)*L^Power_even_trunc(end));
sigma_w=sqrt(sigma_w_squared);

%% State-space matrices

n_state=length(poly_coeff)-1;

Fc=zeros(n_state);
c_coeff=flip(-poly_coeff(2:end));
Fc(end,:)=c_coeff;
Fc(1:end-1,2:end)=eye(n_state-1);

Lc=zeros(n_state,1); Lc(end)=1;
Hc=zeros(1,n_state); Hc(1)=1;

%% Two-sided PSD

S_trunc=[];

if ~isempty(omega_axis)
	Coeff_PSD_even=A_coeff_even_trunc.*L.^Power_even_trunc; 
	Coeff_PSD_all=[Coeff_PSD_even;zeros(size(Coeff_PSD_even))]; Coeff_PSD_all=Coeff_PSD_all(1:end-1);

	S_trunc(1,1,:)=1./polyval(flip(Coeff_PSD_all),omega_axis)*sigma_p^2*L/sqrt(2*pi);
end
