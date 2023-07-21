function [x_k_k,x_k_kmin,x_k_N,P_k_k,P_k_kmin,P_k_N]=KF_RTS(A,B,G,J,Q,R,S,y,p_det,x0,P_0_0,varargin)

%% Kalman Filter and Rauch–Tung–Striebel smoother
%
% Model:
% x(k+1)=A*x(k)+B*p(k)+w(k);
% y(k)=G*x(k)+J*p(k)+v(k);
%
% Inputs:
% A: state matrix
% B: input matrix
% G: output matrix
% J: direct transmission matrix
% Q: state noise covariance
% R: output noise covariance
% S: mixed noise covariance
% y: output vector
% p_det: known (deterministic) input
% x0: initial state estimate
% P_0_0: initial state error covariance
%
% Outputs:
% x_k_k: filter estimate
% x_k_kmin: prediction estimate
% x_k_N: smoothing estimate
% P_k_k: filter error cov
% P_k_kmin: prediction error cov
% P_k_N: smoothing error cov

%% Parse inputs

p=inputParser;
p.KeepUnmatched=true;

addParameter(p,'steadystate',true,@islogical)
addParameter(p,'showtext',true,@islogical)
addParameter(p,'noscaling',false,@islogical)

parse(p,varargin{1:end});

steadystate=p.Results.steadystate;
showtext=p.Results.showtext;
noscaling=p.Results.noscaling;


%%

nx=size(A,1);
ny=size(G,1);
nt=size(y,2);

%%

if isempty(B) & isempty(J)
    B=zeros(nx,1);
    J=zeros(ny,1);
    p=zeros(1,nt);
    
end

%% Transformation in case S is not zero
%
% x(k+1)=A*x(k)+B*p(k)+w(k)
% y(k)=G*x(k)+J*p(k)+v(k)
%
% State equation, add and subtract output term:
% x(k+1)=A*x(k)+B*p(k)+w(k)+S*inv(R)*y(k)-S*inv(R)*y(k)
% x(k+1)=A*x(k)+B*p(k)+w(k)+S*inv(R)*y(k)-S*inv(R)*(G*x(k)+J*p(k)+v(k))
%
% Resulting model:
% x(k+1)=[A-S*inv(R)*G]*x(k)+[B-S*inv(R)*J , S*inv(R)]*[p(k) ; y(k)]+w(k)-S*inv(R)*v(k)
% y(k)=G*x(k)+[J , 0]*[p(k) ; y(k)]+v(k)
%
% The stacked vector [p(k) ; y(k)] is now treated as the deterministic (known) input
% The new process noise is w_star(k)=w(k)+S*inv(R)*v(k)
% The new output noise is v_star(k)=v(k)
% The covariance between w_star(k) and v_star(k) is zero (S_star is zero)
% The states x are not transformed in any way

% State space matrices

if any(any(S))

    A_star=A-S/R*G;
    B_star=[B-S/R*J S/R];
    G_star=G;
    J_star=[J zeros(ny,ny)];
    
    p_star=[p_det ; y];
    Q_star=Q-S/R*S.';
    R_star=R;
    S_star=zeros(size(A,1),size(G,1));
else

    A_star=A;
    B_star=B;
    G_star=G;
    J_star=J;
    
    p_star=p_det;
    Q_star=Q;
    R_star=R;
    S_star=S;
end

%%

[x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KF(A_star,B_star,G_star,J_star,Q_star,R_star,S_star,y,p_star,x0,P_0_0);

[x_k_N,P_k_N]=RTSSmoother(A_star,x_k_k,x_k_kmin,P_k_k,P_k_kmin);


