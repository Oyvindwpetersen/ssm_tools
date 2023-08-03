function [x_filt p_filt Px_k_k_ss Pp_ss M_ss K_ss x_pred Px_k_kmin_ss]=JIS_ss_old(A,B,G,J,y,x0,Q,R,S,P01,varargin)
%% Joint input and state estimation for linear system
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
% M_ss: matrix
% L_ss: matrix
%

%% Parse inputs

p=inputParser;
addParameter(p,'trunc',false,@islogical)
addParameter(p,'scale',false,@islogical)
addParameter(p,'minsteps',100,@isnumeric)
addParameter(p,'maxsteps',100e3,@isnumeric)
addParameter(p,'convtol',1e-6,@isnumeric)
addParameter(p,'dispconv',false,@islogical)
addParameter(p,'Bu',[],@isnumeric)
addParameter(p,'Ju',[],@isnumeric)
addParameter(p,'u',[],@isnumeric)

parse(p,varargin{:});

trunc=p.Results.trunc;
scale=p.Results.scale;
minsteps=p.Results.minsteps;
maxsteps=p.Results.maxsteps;
convtol=p.Results.convtol;
dispconv=p.Results.dispconv;
Bu=p.Results.Bu;
Ju=p.Results.Ju;
u=p.Results.u;

%%

warning('Old version, wrong handling of S');

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

if isempty(x0) | x0==0
    x0=zeros(ns,1);
end

if isempty(P01) | P01==0
    [~,~,~,P01]=KF(A,B,G,J,Q,R,S,y(:,1:min(100,nt)),zeros(np,min(100,nt)),[],[],'noscaling',false,'showtext',false);
end

%assign initial values
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
    Nk=A*Kk*(eye(ny)-J*Mk)+B*Mk;
    P_k_kmin=[A B]*[ P_k_k Pxpkk ; Ppxkk Pp_k_k]*[A.';B.']+Q-Nk*S.'-S*Nk.'; P_k_kmin=forcesym(P_k_kmin);
    
    % Trace
    ratio_trace_Px(k)=trace(abs(P_k_k-P_k_k_prev))./trace(P_k_k);
    ratio_trace_Pp(k)=trace(abs(Pp_k_k-Pp_k_k_prev))./trace(Pp_k_k);

    P_k_k_prev=P_k_k;
    Pp_k_k_prev=Pp_k_k;
    
    if dispconv & mod(k,100)==0
        disp(['***** Step ' num2str(k,'%3.0f') ', ratio_trace_P ' num2str(ratio_trace_Px(k),'%0.3e') ', ratio_trace_Pp ' num2str(ratio_trace_Pp(k),'%0.3e')]);
    end
    
    if k>minsteps & abs(ratio_trace_Px(k)) < convtol & abs(ratio_trace_Pp(k)) < convtol
        convreached=true;
        disp(['Trace convergence reached, k=' num2str(k)]);
        M_ss=Mk;
        K_ss=Kk;
        Px_k_k_ss=P_k_k;
        Pp_ss=Pp_k_k;
        Px_k_kmin_ss=P_k_kmin;
    elseif k>=maxsteps
        figure(); hold on; grid on;
        plot(ratio_trace_Px);
        plot(ratio_trace_Pp);
        set(gca,'YScale','log'); title('Ratio Px');
        error(['Divergence reached, k=' num2str(k)]);
    end
    
end

%% Filter estimates

for k=1:nt
    
    p_filt(:,k)=M_ss*(y(:,k)-G*x_pred(:,k)-Ju*u(:,k));
    x_filt(:,k)=x_pred(:,k)+K_ss*(y(:,k)-G*x_pred(:,k)-J*p_filt(:,k)-Ju*u(:,k));
    x_pred(:,k+1)=A*x_filt(:,k)+B*p_filt(:,k)+Bu*u(:,k);
    
end
x_pred=x_pred(:,1:end-1);

telapsed=toc(tstart);
disp(['JIS calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed*10^3./nt) ' seconds per 1k steps']);

if scale==true
    x_filt=T1*x_filt;
    x_pred=T1*x_pred;
    Px_k_k_ss=T1*Px_k_k_ss*T1.';
    Px_k_kmin_ss=T1*Px_k_kmin_ss*T1.';
    % Set these to empty for safety, not yet checked how these would be affected
    M_ss=[];
    K_ss=[];
end


