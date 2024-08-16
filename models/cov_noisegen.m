function [w,v]=cov_noisegen(Q,R,S,t)

%% Convert covariance from cont to disc 
%
% Model:
% w(k);
% v(k);
%
% Inputs:
% Q: covariance matrix of process noise
% R: covariance matrix of measurement noise
% S: covariance matrix
% t: time vector
%
% Outputs:
% w: process noise
% v: measurement noise

%% Generate noise
nt=length(t);

nq=size(Q,1);
nr=size(R,1);

if isempty(Q)
    Q=zeros(1,1);
    nq=size(Q,1);
end

if isempty(R)
    R=zeros(1,1);
    nr=size(R,1);
end

if isempty(S)
    S=zeros(nq,nr);
end

Ctot=[Q S ; S.' R];

Ctot=(Ctot+Ctot.')/2;

[v,d]=eig(Ctot);
dd=diag(d);

if any(dd<0)

    idx_neg=find(dd<0);
    
    if all( abs(dd(idx_neg))<eps )

        dd(idx_neg)=dd(idx_neg)+eps;

        d_new=diag(dd);

        Ctot=v*d_new*v.';
        
    else
        disp(dd);
    end
end

[wv]=mvnrnd(zeros(size(Ctot,1),1),Ctot,nt).';

range_w=[1:size(Q,1)];
range_v=range_w(end)+[1:size(R,1)];

% Split
w=wv(range_w,:);
v=wv(range_v,:);

