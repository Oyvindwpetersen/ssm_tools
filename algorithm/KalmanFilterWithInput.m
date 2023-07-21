function [x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KalmanFilterWithInput(A,B,G,J,Q,R,S,y,p_det,x0,P_0_0,varargin)

%%

warning('Obsolete, use KF.m instead');

[x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KF(A,B,G,J,Q,R,S,y,p_det,x0,P_0_0,varargin{:});