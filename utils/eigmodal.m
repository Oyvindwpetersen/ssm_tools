function [omega,phi]=eigmodal(K,M,varargin)

%% Find index of x that most closely matches x_val

% Inputs:
% x: vector with all values
% x_val: vector with desired values

% Outputs:
% ind: closest index
% dist: distances/residuals

%%

normalize='mass';

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
        
        alfa=1/sqrt(V(:,i)'*M*V(:,i));
        
        phi(:,i)=V(:,i)*alfa;
        
    end
    
end