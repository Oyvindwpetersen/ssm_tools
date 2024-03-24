function [anti_omega,k0,k1,K]=antires(omega,xi,phi,idx_d,idx_p)

%%

if ~isvector(xi)
    xi=diag(xi);
end

if ~isvector(omega)
    omega=diag(omega);
end

nm=length(omega);

%% Common denominator

coeff_sum=[];

for m=1:nm

    % Start value
    coeff_all=1;

    for n=1:nm

        if n==m
            continue;
        end

        % Polynomial coeffcients for the s^2 polynomial
        coeff_poly2=[1 2*xi(n)*omega(n) omega(n)^2];

        % Polynomial coeffcients the multiplied polynomial (by conv)
        coeff_all=conv(coeff_all,coeff_poly2);

    end

    phi_dm=phi(idx_d,m);
    phi_pm=phi(idx_p,m);

    coeff_sum(m,:)=coeff_all*phi_dm*phi_pm;

end

%% Factorization

% Sum all m=[1:nm]
coeff_sum=sum(coeff_sum,1);

% Normalize such that coefficient of highest term is unity
K=coeff_sum(1);
coeff_sum=coeff_sum/K;

% Roots will be complex pairs
r=roots(coeff_sum);
r=cplxpair(r);

% step=nm-1;

% Select one of each pair
r_sel=r(1:2:end);

for k=1:length(r_sel)

    % Construct polynomial from (s-r)(s-r*)
    p=poly([r_sel(k),conj(r_sel(k))]);

    k1(k)=p(end-1);
    k0(k)=p(end);
end

anti_omega=sqrt(k0);


