function [x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KF(A,B,G,J,Q,R,S,y,p_det,x0,P_0_0,varargin)

%% Kalman filter with known input
%
% Model:
% x(k+1)=A*x(k)+B*p(k)+w(k)
% y(k)=G*x(k)+J*p(k)+v(k)
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

%% Zero matrices

p=p_det;

nx=size(A,1);
ny=size(G,1);
nt=size(y,2);
np=size(p,1);

if isempty(p)
    p=zeros(1,nt);
    np=size(p,1);
end

if isempty(B)
    B=zeros(nx,np);
end

if isempty(J)
    J=zeros(ny,np);
end

if isempty(S)
    S=zeros(nx,ny);
end

if isempty(x0)
    x0=zeros(nx,1);
end

x_hat_k_kmin=zeros(nx,nt);
x_hat_k_k=zeros(nx,nt);
e_k=zeros(ny,nt);

trace_P_k_k=zeros(1,nt);
delta_trace_P_k_k=zeros(1,nt);

if isempty(P_0_0);
    P_0_0=eye(nx);
end

%% Conventional

if steadystate==false
    P_k_kmin=zeros(nx,nx,nt);
    P_k_k=zeros(nx,nx,nt);

    t0=tic;
    for k=1:nt

        %From prev
        if k==1
            x_hat_k_kmin(:,k)=x0;
            P_k_kmin(:,:,k)=P_0_0;
        end

        % Measurement update
        Omega_k=G*P_k_kmin(:,:,k)*G.'+R; Omega_k=forcesym(Omega_k);
        Omega_k_inv=eye(size(Omega_k))/Omega_k;

        P_k_k(:,:,k)=P_k_kmin(:,:,k)-P_k_kmin(:,:,k)*G.'*Omega_k_inv*G*P_k_kmin(:,:,k); P_k_k(:,:,k)=forcesym(P_k_k(:,:,k));

        x_hat_k_k(:,k)=x_hat_k_kmin(:,k)+(P_k_kmin(:,:,k)*G.')*Omega_k_inv*(y(:,k)-G*x_hat_k_kmin(:,k)-J*p(:,k));

    	K_k=(A*P_k_kmin(:,:,k)*G.'+S)*Omega_k_inv;

        % Time update
        %     P_k_kmin(:,:,k+1)=(F-K_k*H)*P_k_kmin(:,:,k)*(F-K_k*H).'+Q+K_k*R*K_k.'-S*K_k.'-K_k*S.'; P_k_kmin(:,:,k+1)=forcesym(P_k_kmin(:,:,k+1));
     	P_k_kmin(:,:,k+1)=A*P_k_kmin(:,:,k)*A.'-(A*P_k_kmin(:,:,k)*G.'+S)*Omega_k_inv*(A*P_k_kmin(:,:,k)*G.'+S).'+Q; P_k_kmin(:,:,k+1)=forcesym(P_k_kmin(:,:,k+1));

        e_k(:,k)=y(:,k)-G*x_hat_k_kmin(:,k)-J*p(:,k);

        x_hat_k_kmin(:,k+1)=A*x_hat_k_kmin(:,k)+B*p(:,k)+K_k*e_k(:,k);

    end

    x_hat_k_kmin(:,end)=[];
    telapsed=toc(t0);

end


%% Steady state

if steadystate==true

    if noscaling==true
        [P_k_kmin_ss,~,~,info]=idare(A.',G.',Q,R,S,'noscaling');
    else
        [P_k_kmin_ss,~,~,info]=idare(A.',G.',Q,R,S);
    end

    if info.Report~=0

        if info.Report==1
            disp('***** DARE solution accuracy poor, running with scaling');
            warning('***** DARE solution accuracy poor, running with scaling');
        end

        if info.Report==2
            warning('***** DARE solution not finite, running with scaling');
        end

    	if info.Report==3
            warning('***** DARE solution not found, running with scaling');
        end

        [P_k_kmin_ss,~,~,info]=idare(A.',G.',Q,R,S);

    end

    if info.Report==1
        disp('***** DARE solution accuracy poor');
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

    % Possibly do further iterations on steady state
    % for k=1:100
    %
    %     % measurement update
    %     Omega_k_ss=G*P_k_kmin_ss*G.'+R; Omega_k_ss=forcesym(Omega_k_ss);
    %     Omega_k_ss_inv=eye(size(Omega_k_ss))/Omega_k_ss;
    %
    %     P_k_k_ss=P_k_kmin_ss-P_k_kmin_ss*G.'*Omega_k_ss_inv*G*P_k_kmin_ss; P_k_k_ss=forcesym(P_k_k_ss);
    % 	K_k_ss=(A*P_k_kmin_ss*G.'+S)*Omega_k_ss_inv;
    %
    %     % time update
    %  	P_k_kmin_ss=A*P_k_kmin_ss*A.'-(A*P_k_kmin_ss*G.'+S)*Omega_k_ss_inv*(A*P_k_kmin_ss*G.'+S).'+Q; P_k_kmin_ss=forcesym(P_k_kmin_ss);
    %
    % end

    t0=tic;
    Mat_precalc=(P_k_kmin_ss*G.')*Omega_k_ss_inv;
    for k=1:nt

        % Initial
        if k==1
            x_hat_k_kmin(:,k)=x0;
        end

        e_k(:,k)=y(:,k)-G*x_hat_k_kmin(:,k)-J*p(:,k);

        % Measurement update
        %     x_hat_k_k(:,k)=x_hat_k_kmin(:,k)+Mat_precalc*(y(:,k)-H*x_hat_k_kmin(:,k));
        x_hat_k_k(:,k)=x_hat_k_kmin(:,k)+Mat_precalc*e_k(:,k); %faster (skipping repeated calculations of e)

        % Time update
        %     x_hat_k_kmin(:,k+1)=F*x_hat_k_kmin(:,k)+K_k_ss*(y(:,k)-H*x_hat_k_kmin(:,k));
        x_hat_k_kmin(:,k+1)=A*x_hat_k_kmin(:,k)+B*p(:,k)+K_k_ss*e_k(:,k); %faster (skipping repeated calculations of e)

    end

    x_hat_k_kmin(:,end)=[];
    telapsed=toc(t0);

end
%%

if showtext==true
    disp(['Kalman filter calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed*10^5./nt) ' seconds per 1M steps']);
end

%% Output

x_k_k=x_hat_k_k;
x_k_kmin=x_hat_k_kmin;

if steadystate==true
    P_k_k=P_k_k_ss;
    P_k_kmin=P_k_kmin_ss;
elseif steadystate==false
    %
end

%% Source

% https://zone.ni.com/reference/en-XX/help/371894J-01/lvsim/sim_disckalmanfilter/
% https://zone.ni.com/reference/en-XX/help/371894J-01/lvcdsimshrd/model_definitions/

% https://se.mathworks.com/help/control/ref/ss.kalman.html
