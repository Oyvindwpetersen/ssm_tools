function [M,K]=mass_spring(m,k)

%% Simple mass spring system
%
% Inputs:
% m: vector with lumped masses
% k: vector with spring stiffness
%
% Outputs:
% M: mass matrix
% K: stiffness matrix
%

%%

M=diag([m]);
K=diag(k+[k(2:end) 0])-diag(k(2:end),-1)-diag(k(2:end),1);
