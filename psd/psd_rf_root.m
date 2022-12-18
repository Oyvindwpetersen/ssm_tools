function [S_rf,rn,n,rd,d]=psd_rf_root(omega,alpha,z_re,z_im,p_re,p_im,forcezero)

%% Spectral density from rational function, with inputs as real and imaginary part of conjugate roots (zeros and poles)
%
% S(omega)=alpha*N(omega)/D(omega)
%
% N(omega)=n0+n2*omega^2+n4*omega^4+...+n_2M*omega^(2M)
% D(omega)=d0+d2*omega^2+d4*omega^4+...+d_2K*omega^(2K)
%
% Inputs:
% omega: frequency vector
% alpha: scale factor
% z_re: real part roots of N-polynomial
% z_im: imaginary part roots of N-polynomial
% p_re: real part roots of D-polynomial
% p_im: imaginary part roots of D-polynomial
% forcezero: true/false, if true then two zero roots are added so that N(omega=0)=0
%
% Outputs:
% S_rf: [1,1,length(omega)] spectral density
% rn: roots of N-polynomial, complex conjugate pairs
% n: polynomial coefficients of N(omega), even terms only, n=[n_2M ... n2 n0]
% rd: roots of D-polynomial, complex conjugate pairs
% d: polynomial coefficients of D(omega), even terms only, d=[d_2K ... d2 d0]
%
%%

if nargin==6
    forcezero=false;
end

%%

% Produce complex conjugate pairs of roots
[rn,n]=pos_poly(z_re,z_im,forcezero);
[rd,d]=pos_poly(p_re,p_im);

% Spectral density
S_rf=psd_rf(omega,n,d,alpha);

% Equivalent but slow, root form of polynomial
% S_rf2(1,1,:)=alpha*prod(omega-rn,1)./ prod(omega-ra,1);

