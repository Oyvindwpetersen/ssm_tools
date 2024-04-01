function [omega_a,xi_a,omega_a0,xi_a0,k0,k1,K]=antires(omega,xi,phi,idx_d,idx_p)

%% Anti resonance frequencies for modal systems
%
% Inputs:
% omega: undamped natural frequencies
% xi: damping ratios
% phi: mode shape matrix
% idx_d: row index of phi for output
% idx_p: row index of phi for input
%
% Outputs:
% omega_a: anti resonance frequencies
% xi_a: anti resonance damping ratios
% omega_a0: anti resonance frequencies (non-conjuage pairs)
% xi_a0: anti resonance damping ratios (non-conjuage pairs)
% k0: factorized coefficients
% k1: factorized coefficients
% K: normalization constant
%

%% Input

if ~isvector(xi)
    xi=diag(xi);
end

if ~isvector(omega)
    omega=diag(omega);
end

nm=length(omega);

%% Common nominator

coeff_sum_tmp=[];

for m=1:nm

    % Start value
    coeff_all=1;

    for n=1:nm

        if n==m
            continue;
        end

        % Polynomial coefficients for the s^2 polynomial
        coeff_poly2=[1 2*xi(n)*omega(n) omega(n)^2];

        % Polynomial coeffcients the multiplied polynomial (by conv)
        coeff_all=conv(coeff_all,coeff_poly2);

    end

    phi_dm=phi(idx_d,m);
    phi_pm=phi(idx_p,m);

    coeff_sum_tmp(m,:)=coeff_all*phi_dm*phi_pm;

end

%% Factorization

% Sum all m=[1:nm]
coeff_sum=sum(coeff_sum_tmp,1);

% If coefficient of highest term is zero, remove it
tol_ph=1e-12;
if (coeff_sum(1))<tol_ph
    coeff_sum=coeff_sum(2:end);
    isreduced=true;
else
    isreduced=false;
end

% Normalize such that coefficient of highest term is unity
K=coeff_sum(1);
coeff_sum=coeff_sum/K;

% Find zeros
lambda_z=roots(coeff_sum);
anti_omega=abs(lambda_z);
anti_xi=-real(lambda_z)./anti_omega;

% Most roots will be complex pairs
r=roots(coeff_sum);
[~,idx_r]=uniquetol(abs(r),1e-6);
r_uni=r(idx_r);

% Find those that correspond to non-conjugate pairs
xi_tmp=anti_xi(idx_r);
bool_nonconj=abs(xi_tmp-1)<1e-6;

omega_a=anti_omega(idx_r(~bool_nonconj));
omega_a0=anti_omega(idx_r(bool_nonconj));

xi_a=anti_xi(idx_r(~bool_nonconj));
xi_a0=anti_xi(idx_r(bool_nonconj));

% If polynomial is not altered, calculate k0 and k1
if isreduced==false

    for k=1:length(r_uni)

        % Construct polynomial from (s-r)(s-r*)
        p=poly([r_uni(k),conj(r_uni(k))]);

        k1(k)=p(end-1);
        k0(k)=p(end);
    end

else
    k1=[];
    k0=[];
end

% anti_omega=sqrt(k0);
% anti_xi=k1./(2*sqrt(k0));

return
%% Debug

clc
close all

H1=0;
omega_axis=linspace(0,40,1e4);
s=1i*omega_axis;
for m=1:nm
    H1=H1+phi(idx_d,m)*phi(idx_p,m)./...
        (s.^2+2*xi(m)*omega(m)*s+omega(m).^2);
end

H1=diag3d(H1);

H1=0;
omega_axis=linspace(0,40,1e4);
s=1i*omega_axis;
for m=1:nm
    H1=H1+phi(idx_d,m)*phi(idx_p,m)./...
        (s.^2+2*xi(m)*omega(m)*s+omega(m).^2);
end

H1=diag3d(H1);

H2_tmp1=1;
for n=1:nm-1
    H2_tmp1=H2_tmp1.*(s.^2+k1(n)*s+k0(n));
end

H2_tmp2=1;
for m=1:nm
    H2_tmp2=H2_tmp2.*(s.^2+2*xi(m)*omega(m)*s+omega(m).^2);
end

H2=K*H2_tmp1./H2_tmp2; %*1e-12;

H2=diag3d(H2)

H3_tmp1=0;
for m=1:nm
    H3_tmp1=H3_tmp1+polyval(coeff_sum_tmp(m,:),s);
end

H3=H3_tmp1./H2_tmp2; %*1e-12;

H3=diag3d(H3)


plottf(omega_axis,abs(H1),abs(H2),abs(H3),'xlim',[0 40],'linestyle',{'-' '--' '-.'});
%%%%


