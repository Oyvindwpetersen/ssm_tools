function ssmod_psd_rf(n,a)

%% Design of state-space model with target PSD
%
% ds/dt=Fc*s(t)+Lc*w(t)
% p(t)=Hc*s(t)
%
% The PSD of the output p(t) is the rational function
%
% S(omega)=N(omega)/A(omega)=alpha*N(omega)/D(omega)
%
% N(omega)=n0+n2*omega^2+n4*omega^4+...
% A(omega)=a0+a2*omega^2+a4*omega^4+a6*omega^6+...
%
% Inputs:
% n: vector with polynomial coefficients n=[n0 n2 n4 ...];
% a: vector with polynomial coefficients a=[a0 a2 a4 a6 ...];
%
% Outputs:
% Fc: state matrix (cont)
% Lc: input matrix (cont)
% Hc: output matrix (cont)
% sigma_w: standard deviation of w
% alpha: scale factor such that d(end)=1

%% Example numbers

% omega=[0:0.01:20];

% a=[1 -3 5 2];
% n=[0 1 0];
% S_exact(1,1,:)=omega.^2./(1-3*omega.^2+5*omega.^4+2*omega.^6);

% a=[0.3^2 1]*2*pi/(2.2^2);
% n=1
% S_exact(1,1,:)=2.2.^2./(2*pi*(0.3.^2+omega.^2));

%%

% Ensure row vector
if size(a,1)>size(a,2); a=a.'; end
if size(n,1)>size(n,2); n=n.'; end

% Scale factor
alpha=1/a(end);

% Highest d coefficient equal to one
d=a*alpha;

% White noise variance
sigma_w_squared=alpha*2*pi;

sigma_w=sqrt(sigma_w_squared);

%% Order of polynomials

num_d=length(d);

% Cut at highest polynomial order not equal to zero
ind_cut=max(find(n~=0));
n=n(1:ind_cut);

num_n=length(n);

if num_n>=num_d
   warning('Order of nominator must be less than denominator');
   warning('E.g. if denominator has omega^6, the max order of nominator is omega^4');
   error('Cannot continue');
end

%% Denominator

% Fill zero for odd terms, reverse every other term
d_all=[];
for k=1:length(d)
    if mod(k,2)==1; signcoeff=1;
    else signcoeff=-1;
    end
    d_all=[d_all signcoeff*d(k) 0];
end
d_all=d_all(1:end-1);

% Find roots
p_all=roots(flip(d_all));

% Identify zero roots, remove temporarily
ind_p_zero=find(p_all==0);
p_all(ind_p_zero)=[];

% Select stable
p_selected=p_all(real(p_all)<0);

% Add back zeros
if ~isempty(ind_p_zero)
    p_zero=zeros(length(ind_p_zero)/2,1);
    p_selected=[p_selected ; p_zero];
end

c_coeff=flip(poly(p_selected)); c_coeff=c_coeff(1:end-1);

%% Nominator

% Fill zero for odd terms, reverse every other term
n_all=[];
for k=1:length(n)
    if mod(k,2)==1; signcoeff=1;
    else signcoeff=-1;
    end
    n_all=[n_all signcoeff*n(k) 0];
end
n_all=n_all(1:end-1);

% Find roots
z_all=roots(flip(n_all));

% Identify zero roots, remove temporarily
ind_z_zero=find(z_all==0);
z_all(ind_z_zero)=[];

% Select stable
z_selected=z_all(real(z_all)<0);

% Add back zeros
if ~isempty(ind_z_zero)
    z_zero=zeros(length(ind_z_zero)/2,1);
    z_selected=[z_selected ; z_zero];
end

b_coeff=flip(poly(z_selected))

% Add zeros to b coeff if the nominator poly is lower than the de
b_zeros=zeros(1,num_d-num_n-1)
b_coeff=[b_coeff b_zeros];


%% State-space matrices

ns=length(c_coeff);

Fc=[zeros(ns-1,1) eye(ns-1) ; -c_coeff];

Lc=zeros(ns,1); Lc(end)=1;

Hc=b_coeff;

%% 
% return
% 
% H=ssmod_tf(Fc,Lc,Hc,zeros(1,1),omega);
% 
% Sw=ones(1,1,length(omega))*sigma_w_squared/(2*pi);
% 
% S_ss=mtimes3(H,Sw,H,'nnh');
% 
% close all
% 
% plotSpectrum(omega,S_exact,S_ss,'LineStyleSet',{'-' '--'},'xlim',[0 10]);
