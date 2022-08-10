function [x y]=ssmod_forward(A,B,G,J,F,x0,p)

%% Simulate (forward solution) of state-space system
%
% Inputs:
% A: state matrix
% B: input matrix
% G: output matrix
% J: direct transmission matrix
% F: matrix for first-order-hold, can be set to empty []
% x0: initial state
% p: input vector
%
% Outputs:
% x: state vector
% y: output vector
%
%% 

ns=size(A,1);
nt=length(p);
nd=size(G,1);
np=size(p,1);

x=zeros(ns,nt);
y=zeros(nd,nt);

if length(x0)==1
    if x0==0;
        x0=zeros(size(A,1),1);
    end
end

x(:,1)=[x0];

if F==0
	F=zeros(size(B));
elseif isempty(F)
	F=zeros(size(B));
end

for k=1:nt-1

x(:,k+1)=A*x(:,k)+B*p(:,k)+F*(p(:,k+1)-p(:,k));

end

if isempty(G) & isempty(J)
    y=[];
else
    y=G*x+J*p; 
    %y for final step, making x and y same length
    y(:,nt)=G*x(:,nt)+J*p(:,nt);
end
