function [x_filt p_filt Px_k_k_ss Pp_k_k_ss M_ss K_ss Kbar_ss]=JIS_ss(A,B,G,J,y,x0,Q,R,S,P01,varargin)
%% Joint input and state estimation for linear systems
%
% Updated with correct handling of mixed covariance S
%
% Model
% x(k+1)=A*x(k)+B*p(k)+w(k);
% y(k)=G*x(k)+J*p(k)+v(k);
%
% Inputs:
% A: state matrix
% B: input matrix
% G: output matrix
% J: direct transmission matrix
% y: output vector
% x0: initial state estimate
% Q: state noise covariance
% R: output noise covariance
% S: mixed noise covariance
% P01: initial state error covariance
%
% Outputs:
% x_filt: filter state estimate
% p_filt: filter input estimate
% P_ss: filter state error covariance
% Pp_ss: filter input error covariance
% M_ss: gain matrix for input
% K_ss: gain matrix for filter estimate of state
% Kbar_ss: gain matrix for prediction estimate of state
%
 
%% Parse inputs

p=inputParser;
addParameter(p,'showtext',true,@islogical)
addParameter(p,'dispconv',true,@islogical)
addParameter(p,'trunc',false,@islogical)
addParameter(p,'scale',false,@islogical)
addParameter(p,'minsteps',10,@isnumeric)
addParameter(p,'maxsteps',100e3,@isnumeric)
addParameter(p,'convtol',1e-6,@isnumeric)
addParameter(p,'Bu',[],@isnumeric)
addParameter(p,'Ju',[],@isnumeric)
addParameter(p,'u',[],@isnumeric)

parse(p,varargin{:});

showtext=p.Results.showtext;
dispconv=p.Results.dispconv;
trunc=p.Results.trunc;
scale=p.Results.scale;
minsteps=p.Results.minsteps;
maxsteps=p.Results.maxsteps;
convtol=p.Results.convtol;
Bu=p.Results.Bu;
Ju=p.Results.Ju;
u=p.Results.u;

%%

% Matrix size
nm=NaN;
ns=size(A,1);
nt=size(y,2);
ny=size(G,1);
np=size(B,2);

% Zero matrices
x_pred=zeros(ns,nt);
x_filt=zeros(ns,nt);
p_filt=zeros(np,nt);

% For state steady calculation
P_k_k_prev=nan(ns);
Pp_k_k_prev=nan(np);

ratio_trace_Px=zeros(1,maxsteps);
ratio_trace_Pp=zeros(1,maxsteps);

if trunc
    lambdaRk=zeros(ny,maxsteps);
    lambdaJRkJ=zeros(np,maxsteps);
    lambdaJPpkkJ=zeros(ny,maxsteps);
end

% Initial state zero 
if isempty(x0) | x0==0
    x0=zeros(ns,1);
end

% Initial covariance from KF 
if isempty(P01) | P01==0
    
    q=1e6*eye(np);
    Q_tmp=Q+B*q*B.';
    R_tmp=R+J*q*J.';
    S_tmp=S+B*q*J.';
    
    [~,~,~,P01]=KF(A,[],G,[],Q_tmp,R_tmp,S_tmp,zeros(ny,10),[],[],[],'noscaling',false,'showtext',false);
end

% Assign initial values
x_pred(:,1)=[x0];
P_k_kmin=P01;

% Scale system
if scale==true
    [A,B,G,J,y,Q,R,S,T1,T2,T1_inv,T2_inv]=ssmod_scaleunitcov(A,B,G,J,y,Q,R,S);
    
    % If deterministic force is present, then scale
    if ~isempty(Bu); Bu=T1_inv*Bu; end
    if ~isempty(Ju); Ju=T2_inv*Ju; end
    
end

% Check if deterministic input are present
if all([isempty(Bu) isempty(Ju) isempty(u)])
	Bu=0;
    Ju=0;
    u=zeros(1,nt);
elseif ~all([isempty(Bu) ~isempty(Ju) ~isempty(u)])
    % With deterministic input
    error('Deterministic input not implemented yet');
else
    size(Bu)
    size(Ju)
    size(u)
    error('All or none of Bu, Ju, and u should be empty');
end

%% Steady state

tstart=tic;
convreached=false; k=0;
while convreached==false

    k=k+1;
    
    %input estimation
    Rk=G*P_k_kmin*G.'+R;  Rk=forcesym(Rk);
    
    % Truncation 1
    if trunc==true & ny>nm
        [V1,lambda_orig1]=eig(Rk);
        [lambda1,I1]=sort(diag(lambda_orig1),1,'descend');
        lambdaRk(:,k)=lambda1;
        V1=V1(:,I1);
        tr1=nm;
        Rk_inv=V1(:,1:tr1)*inv(diag(lambda1(1:tr1)))*V1(:,1:tr1)'; Rk_inv=forcesym(Rk_inv);
        Rk=V1(:,1:tr1)*diag(lambda1(1:tr1))*V1(:,1:tr1)';	Rk=forcesym(Rk);
    else
        Rk_inv=eye(ny)/Rk; Rk_inv=forcesym(Rk_inv);
    end
    
    %trunc 2
    if trunc==true & np>nm
        [V2,lambda_orig2]=eig(J.'*Rk_inv*J);
        [lambda2,I2]=sort(diag(lambda_orig2),1,'descend');
        lambdaJRkJ(:,k)=lambda2;
        V2=V2(:,I2);
        tr2=nm;
        JRkJ_inv= V2(:,1:tr2)*inv(diag(lambda2(1:tr2)))*V2(:,1:tr2)'; JRkJ_inv=forcesym(JRkJ_inv);
    else
        JRkJ_inv=eye(np)/(J.'*Rk_inv*J); JRkJ_inv=forcesym(JRkJ_inv);
    end
    
    Mk=JRkJ_inv*J.'*Rk_inv;
    % p_filt(:,k)=Mk*(y(:,k)-G*x_pred(:,k));
    Pp_k_k=JRkJ_inv;  Pp_k_k=forcesym(Pp_k_k);
    
    % Measurement update
    Kk=P_k_kmin*G.'*Rk_inv;
    % x_filt(:,k)=x_pred(:,k)+Kk*(y(:,k)-G*x_pred(:,k)-J*p_filt(:,k));
    
    %trunc 3
    if trunc==true & ny>nm
        [V3,lambda_orig3]=eig(forcesym(J*Pp_k_k*J.'));
        [lambda3,I3]=sort(diag(lambda_orig3),1,'descend');
        lambdaJPpkkJ(:,k)=lambda3;
        V3=V3(:,I3);
        tr3=min(nm,np);
        JPpkkJ=V3(:,1:tr3)*diag(lambda3(1:tr3))*V3(:,1:tr3)';  JPpkkJ=forcesym(JPpkkJ);
    else
        JPpkkJ=(J*Pp_k_k*J.'); JPpkkJ=forcesym(JPpkkJ);
    end
    
    P_k_k=P_k_kmin - Kk*(Rk-JPpkkJ)*Kk.'; P_k_k=forcesym(P_k_k);
    Pxpkk=-Kk*J*Pp_k_k;
    Ppxkk=Pxpkk';
    
    % Time update
    % x_pred(:,k+1)=A*x_filt(:,k)+B*p_filt(:,k);

    Sbar_k=A*P_k_kmin*G.'+S;
    Kbar_k=Sbar_k/Rk+(B-Sbar_k/Rk*J)*Mk;

    P_k_kmin=A*P_k_kmin*A.'+Q-Sbar_k*Kbar_k.'-Kbar_k*Sbar_k.'+Kbar_k*Rk*Kbar_k.'; P_k_kmin=forcesym(P_k_kmin);
    
    % Trace
    ratio_trace_Px(k)=trace(abs(P_k_k-P_k_k_prev))./trace(P_k_k);
    ratio_trace_Pp(k)=trace(abs(Pp_k_k-Pp_k_k_prev))./trace(Pp_k_k);

    P_k_k_prev=P_k_k;
    Pp_k_k_prev=Pp_k_k;
    
    if dispconv & mod(k,100)==0
        disp(['***** Step ' num2str(k,'%3.0f') ', ratio_trace_P ' num2str(ratio_trace_Px(k),'%0.3e') ', ratio_trace_Pp ' num2str(ratio_trace_Pp(k),'%0.3e')]);
    end
    
    if k>=minsteps & abs(ratio_trace_Px(k)) < convtol & abs(ratio_trace_Pp(k)) < convtol
        convreached=true;
        if dispconv
            disp(['Trace convergence reached, k=' num2str(k)]);
        end
        M_ss=Mk;
        K_ss=Kk;
        Kbar_ss=Kbar_k;
        Rbar_ss=Rk;
        Sbar_ss=Sbar_k;
        Px_k_k_ss=P_k_k;
        Pp_k_k_ss=Pp_k_k;
        Px_k_kmin_ss=P_k_kmin;
    elseif k>=maxsteps
        figure(); hold on; grid on;
        plot(ratio_trace_Px);
        plot(ratio_trace_Pp);
        set(gca,'YScale','log'); title('Ratio Px');
        error(['Divergence reached, k=' num2str(k)]);
    end
    
end

%% Estimate

for k=1:nt
    
    p_filt(:,k)=M_ss*(y(:,k)-G*x_pred(:,k));

    e_k=y(:,k)-G*x_pred(:,k)-J*p_filt(:,k);

    x_filt(:,k)=x_pred(:,k)+K_ss*e_k;

    % x_pred(:,k+1)=A*x_pred(:,k)+B*p_filt(:,k)+Sbar_ss/Rbar_ss*e_k; % This should be the same as below
    x_pred(:,k+1)=A*x_pred(:,k)+Kbar_ss*(y(:,k)-G*x_pred(:,k));
    
end
x_pred=x_pred(:,1:end-1);

telapsed=toc(tstart);

if showtext==true
    disp(['JIS calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed*10^3./nt) ' seconds per 1k steps']);
end

%%

if scale==true
    x_filt=T1*x_filt;
    x_pred=T1*x_pred;
    Px_k_k_ss=T1*Px_k_k_ss*T1.';
    Px_k_kmin_ss=T1*Px_k_kmin_ss*T1.';

    % Set these to empty for safety, not yet checked how these would be affected
    M_ss=[];
    K_ss=[];
    Kbar_ss=[];
    Sbar_ss=[];
end
