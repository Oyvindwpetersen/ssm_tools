function [Fc,Lc,Hc,Qc,sigma_w]=ssmod_periodicdecay(omega0,lambda,sigma_p,kernel,ns)

%%
% 
% Decaying kernel
if strcmpi(kernel,'se')
    [Fq,Lq,Hq,sigma_w]=ssmod_squaredexp(1/lambda,sigma_p,ns,[]);
elseif strcmpi(kernel,'matern')
    p_order=ns-1;
    [Fq,Lq,Hq,sigma_w]=ssmod_matern(lambda,p_order,sigma_p);
else
    error('Kernel not recognized')
end

% Periodic kernel
Fp=[0 -omega0*1 ; omega0*1 0];
Lp=eye(2);
Hp=[1 0];

% Matrices by kronecker product

Fc=kron(Fq,eye(2))+kron(eye(size(Fq)),Fp);
Lc=kron(Lq,Lp);
Hc=kron(Hq,Hp);

qj=1;

Qc=kron(sigma_w^2,qj^2*eye(2));

