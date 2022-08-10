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

Ctot=[Q S ; S.' R];

[wv]=mvnrnd(zeros(size(Ctot,1),1),Ctot,nt).';

range_w=[1:size(Q,1)];
range_v=range_w(end)+[1:size(R,1)];

% Split
w=wv(range_w,:);
v=wv(range_v,:);

