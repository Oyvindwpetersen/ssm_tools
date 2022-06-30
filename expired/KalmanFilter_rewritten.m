function [x_k_k,x_k_kmin,P_k_k,P_k_kmin,e_k,K_k_ss]=KalmanFilter(A,B,C,D,Q,R,S,y,x0,P_0_0,varargin)

%%

% Kalman filter:
% x(k+1)=A*x(k)+w(k);
% y=C*x(k)+v(k);

% Inputs:
% A: state matrix
% C: state influence matrix
% Q=cov(w)
% R=cov(v)
% S=cov(w,v)
% y: measurements
% x0: estimate of state at t=0
% P_0_0: uncertainty estimate for x at t=0

% Outputs:
% x_k_k: state estimate (filter)
% x_k_kmin: state estimate (prediction)
% P_k_k: uncertainty covariance for state (filter)
% P_k_kmin: uncertainty covariance for state (prediction)
% e_k: innovation=y-y_hat
% K_k_ss: Kalman gain

%% Input

p=inputParser;
p.KeepUnmatched=true;

addParameter(p,'steadystate','yes',@ischar)
addParameter(p,'showtext','yes',@ischar)
addParameter(p,'forcescaling','no',@ischar)

parse(p,varargin{1:end});

steadystate=p.Results.steadystate;
showtext=p.Results.showtext;
forcescaling=p.Results.forcescaling;

%% Prepare

% Dimensions
ns=size(A,1);
ny=size(C,1);
nt=size(y,2);

% If no input, set B and D to zero
if isempty(p)
    p=zeros(1,nt);
    B=zeros(ny,1);
    D=zeros(ns,1);
    np=1;
else
    np=size(B,1);
end

% Zero matrices
x_hat_k_kmin=zeros(ns,nt);
x_hat_k_k=zeros(ns,nt);
e_k=zeros(ny,nt);
y_hat=zeros(ny,nt);

trace_P_k_k=zeros(1,nt);
delta_trace_P_k_k=zeros(1,nt);

if isempty(P_0_0); P_0_0=eye(ns); end
if isempty(x0); x0=zeros(size(A,1),1); end


%% Conventional Kalman filter

% From https://zone.ni.com/reference/en-XX/help/371894J-01/lvsim/sim_disckalmanfilter/
% From https://zone.ni.com/reference/en-XX/help/371894J-01/lvcdsimshrd/model_definitions/

% From NI website:

% yhat=C*x_pred+D*p
% x_filt=x_pred+M*(y-yhat)
% Omega=C*P_pred*C'+R
% M=P_pred*C'*inv(Omega)
% P_filt=P_pred-M*C*P_pred

% x_pred_next=A*x_pred+B*p+L*(y-yhat)
% L=(A*P_pred*C'+S)*inv(Omega)
% P_pred_next=A*P_pred*A'+Q-L*(A*P_pred*C'+S)'


%% Conventional filter

if ~strcmpi(steadystate,'yes')
    
P_k_kmin=zeros(ns,ns,nt);
P_k_k=zeros(ns,ns,nt);

t0=tic;
for k=1:nt
    
    % From previous step
    if k==1
        x_hat_k_kmin(:,k)=x0;
        P_k_kmin(:,:,k)=P_0_0;
    end
    
    % Measurement update    
    Omega_k=C*P_k_kmin(:,:,k)*C.'+R; Omega_k=forcesym(Omega_k);
    Omega_k_inv=eye(size(Omega_k))/Omega_k; 

    P_k_k(:,:,k)=P_k_kmin(:,:,k)-P_k_kmin(:,:,k)*C.'*Omega_k_inv*C*P_k_kmin(:,:,k); P_k_k(:,:,k)=forcesym(P_k_k(:,:,k));
    
    y_hat(:,k)=C*x_hat_k_kmin(:,k)+D*p(:,k);
    
    x_hat_k_k(:,k)=x_hat_k_kmin(:,k)+(P_k_kmin(:,:,k)*C.')*Omega_k_inv*(y(:,k)-y_hat(:,k));
    
	K_k=(A*P_k_kmin(:,:,k)*C.'+S)*Omega_k_inv;
    
    % Time update
 	P_k_kmin(:,:,k+1)=A*P_k_kmin(:,:,k)*A.'-(A*P_k_kmin(:,:,k)*C.'+S)*Omega_k_inv*(A*P_k_kmin(:,:,k)*C.'+S).'+Q; P_k_kmin(:,:,k+1)=forcesym(P_k_kmin(:,:,k+1));
      
    x_hat_k_kmin(:,k+1)=A*x_hat_k_kmin(:,k)+B*p(:,k)+K_k*(y(:,k)-y_hat(:,k));
    
    e_k(:,k)=y(:,k)-y_hat(:,k);
    
end

x_hat_k_kmin(:,end)=[];
telapsed=toc(t0);

end


%% Steady state filer

if strcmpi(steadystate,'yes')
    
% Try no scaling first
if strcmpi(forcescaling,'no')
    
    if issparse(F) | issparse(C)
    [P_k_kmin_ss,~,~,info]=idare(full(A).',full(C).',Q,R,S,'noscaling');
    else
    [P_k_kmin_ss,~,~,info]=idare(A.',C.',Q,R,S,'noscaling');
    end
    
    if info.Report~=0

        if info.Report==1
        warning('***** DARE solution accuracy poor, running with scaling');
        end

        if info.Report==2
        warning('***** DARE solution not finite, running with scaling');
        end

        if info.Report==3
        warning('***** DARE solution not found, running with scaling');
        end
        
        forcescaling='yes';
%         [P_k_kmin_ss,~,~,info]=idare(A.',C.',Q,R,S);
    end
end

% If necessary, run with scaling
if strcmpi(forcescaling,'yes') 
	if issparse(F) | issparse(C) | issparse(Q) | issparse(R) | issparse(S)
    [P_k_kmin_ss,~,~,info]=idare(full(A).',full(C).',full(Q),full(R),full(S));
    else
    [P_k_kmin_ss,~,~,info]=idare(A.',C.',Q,R,S);
    end
end

if info.Report==1
    warning('***** DARE solution accuracy poor');    
elseif info.Report==2
    error('DARE solution not finite');
elseif info.Report==3
    error('DARE solution not found, ensure that the system has a steady state solution');
end

% Define matrices for steady state
P_k_kmin_ss=forcesym(P_k_kmin_ss);
Omega_k_ss=(C*P_k_kmin_ss*C.'+R); Omega_k_ss=forcesym(Omega_k_ss);
Omega_k_ss_inv=eye(size(Omega_k_ss))/Omega_k_ss; Omega_k_ss_inv=forcesym(Omega_k_ss_inv);
P_k_k_ss=P_k_kmin_ss-P_k_kmin_ss*C.'*Omega_k_ss_inv*C*P_k_kmin_ss; P_k_k_ss=forcesym(P_k_k_ss);
K_k_ss=(A*P_k_kmin_ss*C.'+S)*Omega_k_ss_inv;

% close all

%Possibly do further iterations on steady state
% doIterationsSteadyState=true;
% doIterationsSteadyState=false;
% if doIterationsSteadyState
% 
% P_k_kmin=zeros(ns,ns,nt);
% P_k_k=zeros(ns,ns,nt);
% 
% P_0_0=P_k_k_ss; 
% P_0_1=P_k_kmin_ss; 
% clear P_k_kmin_ss P_k_k_ss Omega_ss Omega_ss_inv K_k_ss
% 
% for k=1:nt
% 
%     if k==nt & delta_trace_P_k_k(k)>1e-6
%         figure(); plot(delta_trace_P_k_k); ylog();
%         error(['No steady state solution, trace ratio= ' num2str(delta_trace_P_k_k(k),'%0.3e') ]);
%     end
%     
% end
% 
% end

t0=tic;
Mat_precalc=(P_k_kmin_ss*C.')*Omega_k_ss_inv;
for k=1:nt
    
    % Initial
    if k==1
        x_hat_k_kmin(:,k)=x0;
    end
    
    % Innovation
    e_k(:,k)=y(:,k)-C*x_hat_k_kmin(:,k)-D*p(:,k);

    % Measurement update
    x_hat_k_k(:,k)=x_hat_k_kmin(:,k)+Mat_precalc*e_k(:,k); %faster (skipping repeated calculations of e)

    % Time update
    x_hat_k_kmin(:,k+1)=A*x_hat_k_kmin(:,k)+B*p(:,k)+K_k_ss*e_k(:,k); %faster (skipping repeated calculations of e)
    
end
x_hat_k_kmin(:,end)=[];
telapsed=toc(t0);

end
%%

if strcmpi(showtext,'yes')
disp(['Kalman filter calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed*10^5./nt) ' seconds per 1M steps']);
end

%% Output

x_k_k=x_hat_k_k;
x_k_kmin=x_hat_k_kmin;

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