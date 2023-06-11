function S_rf=psd_rf(omega,n,d,alpha)

%% Spectral density from rational function
%
% S(omega)=alpha*N(omega)/D(omega)
%
% N(omega)=n0+n2*omega^2+n4*omega^4+...+n_2M*omega^(2M)
% D(omega)=d0+d2*omega^2+d4*omega^4+...+d_2K*omega^(2K)
%
% Inputs:
% omega: frequency vector
% n: polynomial coefficients of N(omega), even terms only, n=[n_2M ... n2 n0]
% d: polynomial coefficients of D(omega), even terms only, d=[d_2K ... d2 d0]
%
% Outputs:
% S_rf: 1*1*N spectral density
%

%%

if nargin==3
    alpha=1;
end

N=polyvaleven(n,omega);
D=polyvaleven(d,omega);

% S_rf=zeros(1,1,length(omega));
S_rf(1,1,:) = alpha*N./D;
