function H=ssmod_tf(A,B,C,D,omega,dt,varargin)

%% Transfer function for state space model (input-output or input-state)
%
% Inputs:
% A: state matrix (cont or disc)
% B: input matrix (cont or disc)
% C: output matrix (cont or disc)
% D: direct transmission matrix (cont or disc)
% omega: frequency vector in rad/s
% dt: time discretization, set to empty [] if cont
%
% Outputs:
% H: transfer function
%
%%

p=inputParser;
addParameter(p,'type','io',@ischar)
parse(p,varargin{:});
type=p.Results.type;

%%

if nargin==5 | isempty(dt)
    dt=[]; time_type='cont';
else
    time_type='disc';
end

np=size(B,2);
ns=size(A,1);
ny=size(C,1);

if isempty(D)
    D=zeros(ny,np);
end

if strcmpi(type,'io')
    if isempty(D) | isempty(C)
        error('C,D must be given');
    end
end


if size(A,1)<50 | issparse(A)
    fast=false;
else
    fast=true;
end

%%

if strcmpi(type,'io')
    
    H=zeros(ny,np,length(omega));
    
    if time_type=='cont'
        
        if fast==false

            for j=1:length(omega)
                H(:,:,j)=C / (-A+1i*omega(j)*eye(ns)) * B+D;
            end

        else 

            [Psi,Lambda_c]=eig(A); n=size(A,1);
    
            Lc=Psi\B;
            Vc=C*Psi;
            
            lambda_c=diag(Lambda_c);
            for j=1:length(omega)
                val_inv=1./(-lambda_c+1i*omega(j));
                mat_inv=sparse(1:n,1:n,val_inv,n,n);
                H(:,:,j)=Vc*mat_inv*Lc+D;
            end

        end
    
    elseif time_type=='disc'
        
        for j=1:length(omega)
            H(:,:,j)=C / (-A+exp(1i*omega(j)*dt)*eye(ns)) * B+D;
        end
        
    end
    
end

%%

if strcmpi(type,'is')
    
    H=zeros(ns,np,length(omega));
    
    if time_type=='cont'
        
        for j=1:length(omega)
            H(:,:,j)=eye(ns) / (-A+1i*omega(j)*eye(ns)) * B;
        end
        
    elseif time_type=='disc'
        
        for j=1:length(omega)
            H(:,:,j)=eye(ns) / (-A+exp(1i*omega(j)*dt)*eye(ns)) * B;
        end
        
    end
    
end
