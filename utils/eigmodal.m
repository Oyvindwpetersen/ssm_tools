function [omega,phi]=eigmodal(K,M,normalize)

%% Find mode shapes
%
% Inputs:
% K: stiffness matrix
% M: mass matrix
%
% Outputs:
% omega: natural frequencies
% phi: 
%
%%

if nargin<3
    normalize='mass';
end

n=length(K);

[V,D]=eig(K,M,'chol');
omega=diag(D.^0.5);

% Sort by increasing
[~,ind]=sort(omega);
omega=omega(ind);
V=V(:,ind);


if strcmpi(normalize,'unit')
    
    for i=1:n
        denom1=max(V(:,i));
        denom2=min(V(:,i));
        if abs(denom1) > abs(denom2)
            denom=denom1;
        else
            denom=denom2;
        end
        phi(:,i)=V(:,i)/denom;
    end
    
end

if strcmpi(normalize,'mass')
    
    for i=1:n
        alfa=1/sqrt(V(:,i).'*M*V(:,i));
        phi(:,i)=V(:,i)*alfa;
    end
    
end