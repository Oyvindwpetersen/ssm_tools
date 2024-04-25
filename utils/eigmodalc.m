function [omega,xi,phi]=eigmodalc(K,C,M,normalize)

%% Find mode shapes
%
% Inputs:
% K: stiffness matrix
% C: damping matrix
% M: mass matrix
%
% Outputs:
% omega: natural frequencies
% xi: damping ratios
% phi: mode shapes
%
%%

if nargin<4
    normalize='mass';
end

Ac=[zeros(size(K)) eye(size(K)) ; -M\K -M\C];

[lambda,psi,omega,xi]=eigA(Ac);

n=size(K,1);

xi=xi(1:2:end);
omega=omega(1:2:end);
phi=psi(1:n,1:2:end);

if strcmpi(normalize,'mass')

    [mpcw,rotrad,phir,mpc]=mpcweighted(phi);

    if all(mpcw>99.999)
        V=real(phir);
        for i=  1:n
            alfa=1/sqrt(V(:,i).'*M*V(:,i));
            phi(:,i)=V(:,i)*alfa;
        end

    else
        warning('Mode shapes not real, could not be mass normalized');
    end
end


