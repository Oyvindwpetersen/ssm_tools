function [x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KalmanFilter(A,G,Q,R,S,y,x0,P_0_0,varargin)


warning('Obsolete, use KF instead');

%% Kalman filter
%
% Model:
% x(k+1)=A*x(k)+w(k);
% y(k)=G*x(k)+v(k);
%
% Inputs:
% A: state matrix
% G: output matrix
% Q: output noise covariance
% R: state noise covariance
% S: mixed noise covariance
% y: output vector
% x0: initial state estimate
% P_0_0: initial state error covariance
%
% Outputs:
% x_k_k: filter state estimate
% x_k_kmin: prediction state estimate
% P_k_k: filter error covariance
% P_k_kmin: prediction error covariance

%% Parse inputs

p=inputParser;
p.KeepUnmatched=true;

addParameter(p,'steadystate','yes',@ischar)
addParameter(p,'showtext','yes',@ischar)
addParameter(p,'forcescaling','no',@ischar)

parse(p,varargin{1:end});

steadystate=p.Results.steadystate;
showtext=p.Results.showtext;
forcescaling=p.Results.forcescaling;

%% Zero matrices

nx=size(A,1);
ny=size(G,1);
nt=size(y,2);

m_hat_k_kmin=zeros(nx,nt);
m_hat_k_k=zeros(nx,nt);
e_k=zeros(ny,nt);

trace_P_k_k=zeros(1,nt);
delta_trace_P_k_k=zeros(1,nt);

if isempty(P_0_0); P_0_0=eye(nx); end
if isempty(x0); x0=zeros(size(A,1),1); end

%% Conventional Kalman filter

if ~strcmpi(steadystate,'yes')
P_k_kmin=zeros(nx,nx,nt);
P_k_k=zeros(nx,nx,nt);

t0=tic;
for k=1:nt
    
    % From previous step
    if k==1
        m_hat_k_kmin(:,k)=x0;
        P_k_kmin(:,:,k)=P_0_0;
    end
    
    % Measurement update    
    Omega_k=G*P_k_kmin(:,:,k)*G.'+R; Omega_k=forcesym(Omega_k);
    Omega_k_inv=eye(size(Omega_k))/Omega_k; 

    P_k_k(:,:,k)=P_k_kmin(:,:,k)-P_k_kmin(:,:,k)*G.'*Omega_k_inv*G*P_k_kmin(:,:,k); P_k_k(:,:,k)=forcesym(P_k_k(:,:,k));
    
    m_hat_k_k(:,k)=m_hat_k_kmin(:,k)+(P_k_kmin(:,:,k)*G.')*Omega_k_inv*(y(:,k)-G*m_hat_k_kmin(:,k));
    
	K_k=(A*P_k_kmin(:,:,k)*G.'+S)*Omega_k_inv;
    
    % Time update
%     P_k_kmin(:,:,k+1)=(A-K_k*G)*P_k_kmin(:,:,k)*(A-K_k*G).'+Q+K_k*R*K_k.'-S*K_k.'-K_k*S.'; P_k_kmin(:,:,k+1)=forcesym(P_k_kmin(:,:,k+1));
 	P_k_kmin(:,:,k+1)=A*P_k_kmin(:,:,k)*A.'-(A*P_k_kmin(:,:,k)*G.'+S)*Omega_k_inv*(A*P_k_kmin(:,:,k)*G.'+S).'+Q; P_k_kmin(:,:,k+1)=forcesym(P_k_kmin(:,:,k+1));
      
    m_hat_k_kmin(:,k+1)=A*m_hat_k_kmin(:,k)+K_k*(y(:,k)-G*m_hat_k_kmin(:,k));
    
    e_k(:,k)=y(:,k)-G*m_hat_k_kmin(:,k);
    
end

m_hat_k_kmin(:,end)=[];
telapsed=toc(t0);

end


%% Steady state

if strcmpi(steadystate,'yes')
    
% Try no scaling first
if strcmpi(forcescaling,'no')
    
    if issparse(A) | issparse(G)
    [P_k_kmin_ss,~,~,info]=idare(full(A).',full(G).',Q,R,S,'noscaling');
%     P_k_kmin_ss=dare_sparse_newton(A,G,Q,R,S); info.Report=0;
    else
    [P_k_kmin_ss,~,~,info]=idare(A.',G.',Q,R,S,'noscaling');
    end
    
    if info.Report~=0

        if info.Report==1
        warning('***** DARE solution accuracy poor, running with scaling');
    %     error('***** DARE solution accuracy poor, running with scaling');
        end

        if info.Report==2
        warning('***** DARE solution not finite, running with scaling');
        end

        if info.Report==3
        warning('***** DARE solution not found, running with scaling');
        end
        
        forcescaling='yes';
%         [P_k_kmin_ss,~,~,info]=idare(A.',G.',Q,R,S);

    end
end

if strcmpi(forcescaling,'yes') 
	if issparse(A) | issparse(G) | issparse(Q) | issparse(R) | issparse(S)
        [P_k_kmin_ss,~,~,info]=idare(full(A).',full(G).',full(Q),full(R),full(S));
    else
        [P_k_kmin_ss,~,~,info]=idare(A.',G.',Q,R,S);
    end
end

if info.Report==1
    warning('***** DARE solution accuracy poor');    
elseif info.Report==2
    error('DARE solution not finite');
elseif info.Report==3
    error('DARE solution not found');
end

P_k_kmin_ss=forcesym(P_k_kmin_ss);
Omega_k_ss=(G*P_k_kmin_ss*G.'+R); Omega_k_ss=forcesym(Omega_k_ss);
Omega_k_ss_inv=eye(size(Omega_k_ss))/Omega_k_ss; Omega_k_ss_inv=forcesym(Omega_k_ss_inv);
P_k_k_ss=P_k_kmin_ss-P_k_kmin_ss*G.'*Omega_k_ss_inv*G*P_k_kmin_ss; P_k_k_ss=forcesym(P_k_k_ss);
K_k_ss=(A*P_k_kmin_ss*G.'+S)*Omega_k_ss_inv;

% close all

%Possibly do further iterations on steady state
% doIterationsSteadyState=true;
doIterationsSteadyState=false;
if doIterationsSteadyState

P_k_kmin=zeros(nx,nx,nt);
P_k_k=zeros(nx,nx,nt);

P_0_0=P_k_k_ss; 
P_0_1=P_k_kmin_ss; 
clear P_k_kmin_ss P_k_k_ss Omega_ss Omega_ss_inv K_k_ss

for k=1:nt
    
    %From prev 
    if k==1
        P_k_kmin(:,:,k)=P_0_0;
    end
    
    % measurement update    
    Omega_k=G*P_k_kmin(:,:,k)*G.'+R; Omega_k=forcesym(Omega_k);
    Omega_k_inv=eye(size(Omega_k))/Omega_k; 

    P_k_k(:,:,k)=P_k_kmin(:,:,k)-P_k_kmin(:,:,k)*G.'*Omega_k_inv*G*P_k_kmin(:,:,k); P_k_k(:,:,k)=forcesym(P_k_k(:,:,k));
    
% 	K_k=(A*P_k_kmin(:,:,k)*G.'+S)*Omega_k_inv;

    % time update
%     P_k_kmin(:,:,k+1)=(A-K_k*G)*P_k_kmin(:,:,k)*(A-K_k*G).'+Q+K_k*R*K_k.'-S*K_k.'-K_k*S.'; P_k_kmin(:,:,k+1)=forcesym(P_k_kmin(:,:,k+1));
 	P_k_kmin(:,:,k+1)=A*P_k_kmin(:,:,k)*A.'-(A*P_k_kmin(:,:,k)*G.'+S)*Omega_k_inv*(A*P_k_kmin(:,:,k)*G.'+S).'+Q; P_k_kmin(:,:,k+1)=forcesym(P_k_kmin(:,:,k+1));
      
    trace_P_k_k(k)=trace(P_k_k(:,:,k));
    
    if k>1
    delta_trace_P_k_k(k)=abs((trace_P_k_k(k)-trace_P_k_k(k-1))/trace_P_k_k(k));
        
	if delta_trace_P_k_k(k)<1e-6 & k>10
%         disp(['nt=' num2str(k)]);
        P_k_k_ss=P_k_k(:,:,k); P_k_k_ss=forcesym(P_k_k_ss);
        P_k_kmin_ss=P_k_kmin(:,:,k); P_k_kmin_ss=forcesym(P_k_kmin_ss);
        
        P_k_kmin_ss=forcesym(P_k_kmin_ss);
        Omega_k_ss=(G*P_k_kmin_ss*G.'+R); Omega_k_ss=forcesym(Omega_k_ss);
        Omega_k_ss_inv=eye(size(Omega_k_ss))/Omega_k_ss; Omega_k_ss_inv=forcesym(Omega_k_ss_inv);
        P_k_k_ss=P_k_kmin_ss-P_k_kmin_ss*G.'*Omega_k_ss_inv*G*P_k_kmin_ss; P_k_k_ss=forcesym(P_k_k_ss);
        K_k_ss=(A*P_k_kmin_ss*G.'+S)*Omega_k_ss_inv;

        figure(); 
        plot(delta_trace_P_k_k); ylog;
        break; 
    end
    
    end

    if k==nt & delta_trace_P_k_k(k)>1e-6
        figure(); plot(delta_trace_P_k_k); ylog();
        error(['No steady state solution, trace ratio= ' num2str(delta_trace_P_k_k(k),'%0.3e') ]);
    end
    
end

end

t0=tic;
Mat_precalc=(P_k_kmin_ss*G.')*Omega_k_ss_inv;
for k=1:nt
    
    % Initial
    if k==1
        m_hat_k_kmin(:,k)=x0;
    end
    
    %
    e_k(:,k)=y(:,k)-G*m_hat_k_kmin(:,k);

    % measurement update
%     m_hat_k_k(:,k)=m_hat_k_kmin(:,k)+Mat_precalc*(y(:,k)-G*m_hat_k_kmin(:,k));
    m_hat_k_k(:,k)=m_hat_k_kmin(:,k)+Mat_precalc*e_k(:,k); %faster (skipping repeated calculations of e)

    %time
%     m_hat_k_kmin(:,k+1)=A*m_hat_k_kmin(:,k)+K_k_ss*(y(:,k)-G*m_hat_k_kmin(:,k));
    m_hat_k_kmin(:,k+1)=A*m_hat_k_kmin(:,k)+K_k_ss*e_k(:,k); %faster (skipping repeated calculations of e)
    
end
m_hat_k_kmin(:,end)=[];
telapsed=toc(t0);

end
%%

if strcmpi(showtext,'yes')
disp(['Kalman filter calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed*10^5./nt) ' seconds per 1M steps']);
end

%% Output

x_k_k=m_hat_k_k;
x_k_kmin=m_hat_k_kmin;

if strcmpi(steadystate,'yes')
P_k_k=P_k_k_ss;
P_k_kmin=P_k_kmin_ss;
elseif ~strcmpi(steadystate,'yes')
% P_k_k=P_k_k;
% P_k_kmin_ss=P_k_kmin_ss;
end

end

%%

function B=forcesym(A)

B=(A+A.')/2;

end
