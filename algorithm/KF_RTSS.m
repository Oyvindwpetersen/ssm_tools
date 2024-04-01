function [x_k_k,x_k_kmin,x_k_N,P_k_k,P_k_kmin,P_k_N]=KF_RTSS(A,B,G,J,Q,R,S,y,p_det,varargin)

%% Kalman filter and Rauch–Tung–Striebel smoother
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
%

%% Parse inputs

p=inputParser;
p.KeepUnmatched=true;

addParameter(p,'showtext',true,@islogical)
addParameter(p,'noscaling',false,@islogical)
addParameter(p,'method','standard',@ischar) %standard, decorr
addParameter(p,'x0',[],@isnumerical)
addParameter(p,'P0',[],@isnumerical)

parse(p,varargin{1:end});

showtext=p.Results.showtext;
noscaling=p.Results.noscaling;
method=p.Results.method;
x0=p.Results.x0;
P0=p.Results.P0;

%% Parameters

nx=size(A,1);
ny=size(G,1);
nt=size(y,2);

%% Transformation in case S is not zero
%
% Model:
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
% The new process noise is w_star(k)=w(k)-S*inv(R)*v(k)
% The new output noise is v_star(k)=v(k)
% The covariance between w_star(k) and v_star(k) is zero (S_star is zero)
% The states x are not transformed in any way
%

%% No correlation between process and measurement noise

if ~any(any(S))

    [x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KF(A,B,G,J,Q,R,S,y,p_det,x0,P0,'showtext',showtext,'noscaling',noscaling);

    [x_k_N,P_k_N]=RTSS(A,x_k_k,x_k_kmin,P_k_k,P_k_kmin,'showtext',showtext);

end

%% Correlated noise

if any(any(S))

    if strcmpi(method,'standard')

        % Run KF as standard
        [x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KF(A,B,G,J,Q,R,S,y,p_det,x0,P0,'showtext',showtext,'noscaling',noscaling);

        % Decorrelate for RTS smoother, this only affects the state matrix

        % Optimal and Robust estimation, p 140
        % It is worth remarking that the backward recursive smoother depends neither
        % on the data nor on the deterministic input

        [x_k_N,P_k_N]=RTSS(A-S/R*G,x_k_k,x_k_kmin,P_k_k,P_k_kmin,'showtext',showtext);

    elseif strcmpi(method,'decorr')

        % Decorrelate both for KF and RTS smoother

        A_star=A-S/R*G;
        G_star=G;

        if isempty(B) & isempty(J)
            B_star=[S/R];
            J_star=[zeros(ny,ny)];
        else
            B_star=[B-S/R*J S/R];
            J_star=[J zeros(ny,ny)];
        end

        p_star=[p_det ; y];
        Q_star=Q-S/R*S.';
        R_star=R;
        S_star=zeros(nx,ny);

        [x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KF(A_star,B_star,G_star,J_star,Q_star,R_star,S_star,y,p_star,x0,P0,'showtext',showtext,'noscaling',noscaling);

        [x_k_N,P_k_N]=RTSS(A_star,x_k_k,x_k_kmin,P_k_k,P_k_kmin,'showtext',showtext);

    end

end

