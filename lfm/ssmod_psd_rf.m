function [Fc,Lc,Hc,sigma_w]=ssmod_psd_rf(n,d,alpha)

%% Design of state-space model with target spectral density
%
% Model:
% ds/dt=Fc*s(t)+Lc*w(t)
% p(t)=Hc*s(t)
%
% The (two-sided) PSD of the output p(t) is the rational function
%
% S(omega)=alpha*N(omega)/D(omega)
%
% N(omega)=n0+n2*omega^2+n4*omega^4+...+n_2M*omega^(2M)
% D(omega)=d0+d2*omega^2+d4*omega^4+...+d_2K*omega^(2K)
%
% Inputs:
% n: polynomial coefficients of N(omega), even terms only, n=[n_2M ... n2 n0]
% d: polynomial coefficients of D(omega), even terms only, d=[d_2K ... d2 d0]
%
% Outputs:
% Fc: state matrix (cont)
% Lc: input matrix (cont)
% Hc: output matrix (cont)
% sigma_w: standard deviation of w
%

%% Example numbers
%
% omega=[0:0.01:20];
%
% d=[1 5 -3 8];
% n=[1 0];
% alpha=5;
% S_exact(1,1,:)=5*omega.^2./(8-3*omega.^2+5*omega.^4+omega.^6);
%
%%

% Ensure row vector
if size(d,1)>size(d,2); d=d.'; end
if size(n,1)>size(n,2); n=n.'; end

% Cut at highest polynomial order not equal to zero
ind_first_nonzero=find(n~=0,1);

if ind_first_nonzero>1
    n=n(ind_first_nonzero:end);
end

if abs(n(1)-1)>1e-12
    n
    error('First (i.e. highest order) coefficient in n must be 1');
end

if abs(d(1)-1)>1e-12
    d
    error('First (i.e. highest order) coefficient in d must be 1');
end

if any(isinf(abs(n)))
    n
    error('n coefficients contain inf');
end

if any(isinf(abs(d)))
    d
    error('d coefficients contain inf');
end

% White noise variance
sigma_w_squared=alpha*2*pi;

sigma_w=sqrt(sigma_w_squared);

%% Order of polynomials

num_d=length(d);

if d(1)==0
    error('Highest d coefficient is zero. This means the state-space model is smaller than intended.');
end

num_n=length(n);

if num_n>=num_d
   warning('Order of nominator must be less than denominator');
   warning('E.g. if denominator has omega^6, the max order of nominator is omega^4');
   error('Cannot continue');
end

%% Denominator

% Fill zero for odd terms, flip sign of every other term
signflip=[];
for k=0:(length(d)-1)
    signflip=[(-1)^k signflip] ;
end
d_signflip=d.*signflip;
d_signflip_all=reshape([d_signflip ; zeros(size(d_signflip))],1,[]); d_signflip_all(end)=[];

% Find roots
p_all=roots(d_signflip_all);

% Identify zero roots, remove temporarily
ind_p_zero=find(p_all==0);
p_all(ind_p_zero)=[];

% Select stable
p_selected=p_all(real(p_all)<0);

% Substitute back zero roots
if ~isempty(ind_p_zero)
    p_zero=zeros(length(ind_p_zero)/2,1);
    p_selected=[p_selected ; p_zero];
end

% a_coeff=[a(k),...,a2,a1,a0]
a_coeff=poly(p_selected); 

if abs(a_coeff(1)-1)>1e-6
    a_coeff
    warning('Highest a coefficient should be 1. Something is wrong, check this line');
end

% a_coeff=[a(k-1),...,a2,a1,a0]
a_coeff=a_coeff(2:end); 

%% Nominator

% Fill zero for odd terms, flip sign of every other term
signflip=[];
for k=0:(length(n)-1)
    signflip=[(-1)^k signflip];
end
n_signflip=n.*signflip;
n_signflip_all=reshape([n_signflip ; zeros(size(n_signflip))],1,[]); n_signflip_all(end)=[];

% Find roots
z_all=roots(n_signflip_all);

% Identify zero roots, remove temporarily
ind_z_zero=find(z_all==0);
z_all(ind_z_zero)=[];

% Select stable
z_selected=z_all(real(z_all)<0);

% Substitute back zero roots
if ~isempty(ind_z_zero)
    z_zero=zeros(length(ind_z_zero)/2,1);
    z_selected=[z_selected ; z_zero];
end

% c_coeff=[c(k),...,c2,c1,c0]
c_coeff=poly(z_selected);

% Add zeros to c coeff if the nominator poly is lower than the denominator minus one (e.g. D~omega^6 and N~omega^2)
c_zeros=zeros(1,num_d-num_n-1);

%% State-space model

ns=length(a_coeff);

Fc=[zeros(ns-1,1) eye(ns-1) ; -flip(a_coeff)];

Lc=zeros(ns,1); Lc(end)=1;

Hc=[flip(c_coeff) c_zeros];

