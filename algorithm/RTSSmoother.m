function [x_k_N,P_k_N,P_klag_N]=RTSSmoother(A,x_k_k,x_k_kmin,P_k_k,P_k_kmin,varargin)

%% RTS smoother
%
% Inputs:
% A: state matrix
% x_k_k: filter state estimate from Kalman filter
% x_k_kmin: prediction state estimate from Kalman filter
% P_k_k: filter error covariance from Kalman filter
% P_k_kmin: prediction error covariance from Kalman filter
%
% Outputs:
% x_k_N: smoothed state estimate
% P_k_N: smoothed state error covariance
% P_klag_N: error covariance between k and k-1 given all N data
%
% Note: 
% For S~=0, the system must be transformed so that the process and maesurement noise is uncorrelated. A_star=A-S/R*G;
% See Niu (2011) and Verify_KF_RTS.m 

%% Parse inputs

p=inputParser;

addParameter(p,'steadystate','yes',@ischar)
addParameter(p,'showtext','yes',@ischar)
addParameter(p,'skipRTS','no',@ischar)
addParameter(p,'G',[],@isnumeric)
addParameter(p,'R',[],@isnumeric)

parse(p,varargin{1:end});

steadystate=p.Results.steadystate;
showtext=p.Results.showtext;
skipRTS=p.Results.skipRTS;
G=p.Results.G;
R=p.Results.R;

%% Note on RTS

% From Optimal and Robust estimation, p 140
% It is worth remarking that the backward recursive smoother depends neither
% on the data nor on the deterministic input

%% Check if skip

if strcmpi(skipRTS,'yes')
    x_k_N=NaN*ones(size(x_k_k));
    P_k_N=NaN*ones(size(P_k_k));
    P_klag_N_ss=NaN;
    return
end

if isempty(G)
    skip_Pklag=true;
else
    skip_Pklag=false;
end

%% Zero matrices

m_hat_k_k=x_k_k;
m_hat_k_kmin=x_k_kmin;

nx=size(A,1);
% ny=size(H,1);
nt=size(x_k_k,2);
m_hat_k_N=zeros(nx,nt);
m_hat_k_N(:,nt)=m_hat_k_k(:,nt);

%% Conventional

if ~strcmpi(steadystate,'yes')

    P_k_N=zeros(nx,nx,nt);
    P_k_N(:,:,nt)=P_k_k(:,:,nt);

    if ~skip_Pklag
        Omega=G*P_k_kmin(:,:,nt)*G.'+R;
        P_klag_N=zeros(nx,nx,nt);
        P_klag_N(:,:,nt)=...
        (eye(nx)-P_k_kmin(:,:,nt)*G.'/Omega*G)*A*P_k_k(:,:,nt-1);
    end

    t0=tic;
    for k=(nt-1):-1:1

        N_k=P_k_k(:,:,k)*A.'/P_k_kmin(:,:,k+1);
        m_hat_k_N(:,k)=m_hat_k_k(:,k)+N_k*(m_hat_k_N(:,k+1)-m_hat_k_kmin(:,k+1));
        P_k_N(:,:,k)=P_k_k(:,:,k)+N_k*(P_k_N(:,:,k+1)-P_k_kmin(:,:,k+1))*N_k.';

        if ~skip_Pklag & k>1
        N_kmin=P_k_k(:,:,k-1)*A.'/P_k_kmin(:,:,k);
        P_klag_N(:,:,k)=P_k_k(:,:,k)*N_kmin+N_k*(P_klag_N(:,:,k+1)-A*P_k_k(:,:,k))*N_kmin.';
        end
    end
    telapsed=toc(t0);

end

%% Steady state

if strcmpi(steadystate,'yes')

    P_k_k_ss=P_k_k;
    P_k_kmin_ss=P_k_kmin;

    N_k_ss=P_k_k_ss*A.'/P_k_kmin_ss;

    Q_temp=P_k_k_ss-N_k_ss*P_k_kmin_ss*N_k_ss.'; Q_temp=forcesym(Q_temp);
    P_k_N_ss=dlyap(N_k_ss,Q_temp); P_k_N_ss=forcesym(P_k_N_ss);

    t0=tic;
    for k=(nt-1):-1:1
        m_hat_k_N(:,k)=m_hat_k_k(:,k)+N_k_ss*(m_hat_k_N(:,k+1)-m_hat_k_kmin(:,k+1));
    end
    telapsed=toc(t0);

    P_klag_N_ss=[];
end

%% Calculate the covariance

% This the covariance of the state error between time step k and k-1 given all N time data (smoothing)

% From Cara p. 124. Equation for the Lag-One Covariance Smoother

% At steady state, this covariance should be:
% P_klag_N_ss=P_k_k_ss*N_k_ss.'+N_k_ss(P_klag_N_ss-A*P_k_k_ss)*N_k_ss.';
% Reformulated:
% 0=-P_klag_N_ss    +    N_k_ss*P_klag_N_ss*N_k_ss.'    +    P_k_k_ss*N_k_ss.'    -    N_k_ss*A*P_k_k_ss*N_k_ss.'

% This is a Sylvester equation AXB-X+C=0;
% where
% X=P_klag_N_ss
% A=N_k_ss
% B=N_k_ss.'
% C=P_k_k_ss*N_k_ss.'    -    N_k_ss*A*P_k_k_ss*N_k_ss.'

if strcmpi(steadystate,'yes') & ~skip_Pklag

    A_temp=N_k_ss;
    B_temp=N_k_ss.';
    C_temp=P_k_k_ss*N_k_ss.'    -    N_k_ss*A*P_k_k_ss*N_k_ss.';

    % Solve Sylvester equation with matlab function for Lyapunov
    P_klag_N_ss=dlyap(A_temp,B_temp,C_temp);

else
    

end

%%

if strcmpi(showtext,'yes')
disp(['RTS smoother calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed*10^6./nt) ' seconds per 1M steps']);
end


%% Output

x_k_N=m_hat_k_N;

if strcmpi(steadystate,'yes')
    P_k_N=P_k_N_ss;
    P_klag_N=P_klag_N_ss;
else
    % P_k_N=P_k_N;
    % P_klag_N=P_klag_N;
end


