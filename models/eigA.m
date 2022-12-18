function [lambda,phi,omega,xi]=eigA(A,dt)

%% Eigenvalues for A-matrix
%
% Inputs:
% A: state matrix (cont or disc)
% dt: time discretization, if omitted or empty then A=Ac (continuous time); if provided, then A=Ad
%
% Outputs:
% lambda: eigenvalues of A
% phi: eigenvectors of A
% omega: undamped natural frequency, abs(lambda)
% xi: damping ratio, -real(lambda)./abs(lambda)
%
%% If A is a cell with many A-matrices, calculate for all

% Stack 3D

if iscell(A)
    A_cell=A;  A=zeros(size(A_cell{1},1),size(A_cell{1},2),length(A_cell));
    for k=1:length(A_cell)
        A(:,:,k)=A_cell{k};
    end
end

%% Case single A

if size(A,3)==1

[v,d_discrete]=eig(A);

if nargin==1 | isempty(dt)
% Ac=A;
d_cont=d_discrete;
elseif nargin==2
d_cont=log(d_discrete)./dt;
end

[~,i_sort]=sort(abs(diag(d_cont)));
d_cont=diag(d_cont); d_cont=d_cont(i_sort);
v=v(:,i_sort);

v=v./max(abs(v),[],1);

lambda=d_cont;
phi=v;
omega=abs(lambda); 
xi=-real(lambda)./abs(lambda); 

end

%% Case multiple A

A_all=A;
s3=size(A_all,3);
s1=size(A_all,1);

if s3>1
    
    for k=1:s3
	A=A_all(:,:,k);
        
    [v,d_discrete]=eig(A);

    if nargin==1 | isempty(dt)
    % Ac=A;
    d_cont=d_discrete;
    elseif nargin==2
    d_cont=log(d_discrete)./dt;
    end

    [~,i_sort]=sort(abs(diag(d_cont)));
    d_cont=diag(d_cont); d_cont=d_cont(i_sort);
    v=v(:,i_sort);

    v=v./max(abs(v),[],1);

    lambda=d_cont;
%     phi=v;
    omega=abs(lambda);
    xi=-real(lambda)./abs(lambda);
    
    v_all(:,:,k)=v;
    lambda_all(:,k)=lambda;
    omega_all(:,k)=omega;
    xi_all(:,k)=xi;
    
    end
    
%     v=v_all;
    lambda=lambda_all;
    omega=omega_all;
    xi=xi_all;
    phi=v_all;
end