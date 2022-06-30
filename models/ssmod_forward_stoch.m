function [x y]=ssmod_forward_stoch(A,B,G,J,F,x0,p,w,v)

%% Simulate (forward solution) of state-space system with stochastic noise

% Inputs:
% A: state matrix
% B: input matrix
% G: output matrix
% J: direct transmission matrix
% F: matrix for first-order-hold, can be set to empty []
% x0: initial state
% p: input vector, can be set to empty [] or zero
% w: state noise
% v: output vector

% Outputs:
% x: state vector
% y: output vector

%% 

ns=size(A,1);
nd=size(G,1);
np=size(p,1);
nt=size(w,2);

x=zeros(ns,nt);
y=zeros(nd,nt);

if length(x0)==1 & x0==0
	x0=zeros(size(A,1),1);
end

x(:,1)=[x0];

if F==0 | isempty(F)
	F=zeros(size(B));
end

if all(all(F==0))

	% Case no force
    if all(all(p==0)) | isempty(p)

        for k=1:nt-1
            x(:,k+1)=A*x(:,k)+w(:,k);
            y(:,k)=G*x(:,k)+v(:,k);
        end

        % Final time step
        y(:,nt)=G*x(:,nt)+v(:,nt);

	% Case with force
    else

        for k=1:nt-1
            x(:,k+1)=A*x(:,k)+B*p(:,k)+w(:,k);
            y(:,k)=G*x(:,k)+J*p(:,k)+v(:,k);
        end

        % Final time step
        y(:,nt)=G*x(:,nt)+J*p(:,nt)+v(:,nt);

    end
else
    for k=1:nt-1
        x(:,k+1)=A*x(:,k)+B*p(:,k)+F*(p(:,k+1)-p(:,k))+w(:,k);
        y(:,k)=G*x(:,k)+J*p(:,k)+v(:,k);
    end

    % Final time step
    y(:,nt)=G*x(:,nt)+J*p(:,nt)+v(:,nt);
end


