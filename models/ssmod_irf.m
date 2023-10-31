function [tau,h]=ssmod_irf(A,B,C,D,dt,tau_max,tau_vector)

%% Impulse response function for state space model
%
% Inputs:
% A: state matrix (cont or disc)
% B: input matrix (cont or disc)
% C: output matrix (cont or disc)
% D: direct transmission matrix (cont or disc)
% dt: time discretization, set to empty [] if cont
% tau_max: max time lag in [s] (only for disc)
% tau_vector: time lag vector (only for cont)
%
% Outputs:
% tau: time lag vector
% h: impulse response
%

%%

if nargin==6 
    tau_vector=[];
end

%%

ns=size(A,1);
np=size(B,2);
ny=size(C,1);

if isempty(D)
    D=zeros(ny,np);
end
    
%% Disc time

if ~isempty(dt)

nlag_max=ceil(tau_max/dt);
ind=0;
h=zeros(size(D,1),size(D,2),nlag_max);

tau=[0:nlag_max]*dt;

for k=0:nlag_max
    
    ind=ind+1;
    
    if k==0
        h(:,:,ind)=D;
    else
%         H(:,:,ind)=Cd*Ad^(k-1)*Bd;
        
        if k==1
            temp_CA=C;
        else
            temp_CA=temp_CA*A;
        end
        
        h(:,:,ind)=temp_CA*B;
    end
        
end

end


%% Cont time

if isempty(dt)
    
nlag_max=length(tau_vector);
h=zeros(size(D,1),size(D,2),nlag_max);

% tau=[0:nlag_max]*dt;

for ind=1:nlag_max
    
    h(:,:,ind)=C*expmnorm(A*tau_vector(ind))*B;
    if ind==0
        h(:,:,ind)=h(:,:,ind)+D;
    end
        
end

tau=tau_vector;
end


