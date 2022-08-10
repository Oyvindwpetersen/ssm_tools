function [Fc,Lc,Hc,sigma_w,S_trunc]=ssmod_squaredexp(L,sigma_p,N_terms,w_axis)

%% State-space model approximation of SE kernel

% Inputs:
% L: length scale
% sigma_p: magnitude factor
% N_terms: number of terms in Taylor series, the state-space model will have N_terms-1 states 
% w_axis: axis for spectral density

% Outputs:
% Fc: state kernel matrix in state-space model
% Lc: input matrix in state-space model
% Hc: output matrix in state-space model
% sigma_w: standard deviation of white noise input w(t) in state-space model 
% S_trunc: PSD of approximation model

% CF for SE model:
% Kappa(t)=sigma_p^2*exp(-0.5*t^2/L^2)

% By the inverse Fourier transform:
% S(w)=1/(2*pi) int Kappa(t) exp(-i*w*tau) dtau
% S(w)=exp(-0.5*L^2*w^2)*sigma_p^2*L/sqrt(2*pi)

%% Checks

if N_terms>7 | N_terms<2
    error('N_terms cannot be larger than 7 or smaller than 2');
end

%%

% Taylor expansion of exp(-0.5*L^2*w^2): 1/(w^0*L^0*A0+w^2*L^2*A2+w^4*L^4*A4+...until...w^10*L^10*A10)

% Powers and constants A for the polynomial
Power_even=[0 2 4 6 8 10];
A_coeff_even=[1 1/2 1/8 1/48 1/384 1/3840 1/46080 1/645120];

% Truncate at N terms
Power_even_trunc=Power_even(1:N_terms);
A_coeff_even_trunc=A_coeff_even(1:N_terms);

% Rescale coeff so that highest poly has coeff unity
Coeff_even_trunc_scaled=A_coeff_even_trunc/(A_coeff_even_trunc(end)*L^Power_even_trunc(end));

% Constant value on top of fraction
Constant0=sigma_p^2*L/sqrt(2*pi)/(A_coeff_even_trunc(end)*L^Power_even_trunc(end));

% Roots of poly
% Fill zero coeff for odd w^1,w^3,...
Coeff_all_trunc_scaled=[Coeff_even_trunc_scaled;zeros(size(Coeff_even_trunc_scaled))];
Coeff_all_trunc_scaled=Coeff_all_trunc_scaled(1:end-1);

roots_w=roots(flip(Coeff_all_trunc_scaled));
roots_iw=1i*roots_w;

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
a_coeff=flip(-poly_coeff(2:end))./L.^([n_state:-1:1]);
Fc(end,:)=a_coeff;
Fc(1:end-1,2:end)=eye(n_state-1);

Lc=zeros(n_state,1); Lc(end)=1;
Hc=zeros(1,n_state); Hc(1)=1;

%% Two-sided PSD

S_trunc=[];

if ~isempty(w_axis)
Coeff_PSD_even=A_coeff_even_trunc.*L.^Power_even_trunc; 
Coeff_PSD_all=[Coeff_PSD_even;zeros(size(Coeff_PSD_even))]; Coeff_PSD_all=Coeff_PSD_all(1:end-1);

S_trunc(1,1,:)=1./polyval(flip(Coeff_PSD_all),w_axis)*sigma_p^2*L/sqrt(2*pi);
end
