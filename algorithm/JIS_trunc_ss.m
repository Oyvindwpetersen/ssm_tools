function [x_filt p_filt P_ss Pp_ss M_ss L_ss] = JIS_trunc_ss(A,B,G,J,y,x0,R,Q,S,P01,varargin)

%% Joint input and state estimation for linear system

% Model
% x(k+1)=A*x(k)+B*p(k)+w(k);
% y=G*x(k)+J*p(k)+v(k);

% Inputs:
% A: state matrix
% B: input matrix
% G: output matrix
% J: direct transmission matrix
% y: output vector
% x0: initial state estimate
% R: output noise covariance
% Q: state noise covariance
% S: mixed noise covariance
% P01: initial state error covariance

% Outputs:
% x_filt: filter state estimate
% p_filt: filter input estimate
% P_ss: filter state error covariance
% Pp_ss: filter input error covariance
% M_ss: matrix
% L_ss: matrix


%%
p=inputParser;
addParameter(p,'trunc','no',@ischar)
addParameter(p,'plot','no',@ischar)
addParameter(p,'scale','no',@ischar)
addParameter(p,'minsteps',100,@isnumeric)
addParameter(p,'maxsteps',100e3,@isnumeric)
addParameter(p,'convtol',1e-6,@isnumeric)
addParameter(p,'dispconv',false,@islogical)

parse(p,varargin{:});
doTruncation = p.Results.trunc;
doPlot = p.Results.plot;
doScale = p.Results.scale;
minsteps = p.Results.minsteps;
maxsteps = p.Results.maxsteps;
convtol = p.Results.convtol;
dispconv = p.Results.dispconv;

%%

ns=size(A,1);
nt=length(y);
ny=size(G,1);
np=size(B,2);
nm=NaN;

%zero matrices
x_pred=zeros(ns,nt);
x_filt=zeros(ns,nt);
p_filt=zeros(np,nt);

lambdaRk=zeros(ny,maxsteps);
lambdaJRkJ=zeros(np,maxsteps);
lambdaJPpkkJ=zeros(ny,maxsteps);

%assign initial values
x_pred(:,1)=[x0];
Pkk_=P01;

if strcmpi(doScale,'yes')
%    [A,B,G,J,y,Q,R,S,T1,T2]=scaleStateSpaceModelUnitCov(A,B,G,J,y,Q,R,S); 
   [A,B,G,J,y,Q,R,S,T1,T2,T1_inv,T2_inv]=ssmod_scaleunitcov(A,B,G,J,y,Q,R,S);
end

%% Steady state

clear r_traceP r_tracePp

Pkk_prev=zeros(ns); Pkk_current=zeros(ns);
Ppkk_prev=zeros(np); Ppkk_current=zeros(np);

tstart=tic;
convreached=false; k=0;
while convreached==false    
    k=k+1;
    
    %input estimation
    Rk=G*Pkk_*G'+R;  Rk=forcesym(Rk); 
    
    % Truncation 1
    if 0%ny>nm & strcmpi(doTruncation,'yes');
        [V1,lambda_orig1] = eig(Rk);
        [lambda1,I1] = sort(diag(lambda_orig1),1,'descend');
        lambdaRk(:,k) = lambda1;
        V1 = V1(:,I1);
        tr1 = nm;
        Rk_inv = V1(:,1:tr1)*inv(diag(lambda1(1:tr1)))*V1(:,1:tr1)'; Rk_inv=forcesym(Rk_inv);
        Rk = V1(:,1:tr1)*diag(lambda1(1:tr1))*V1(:,1:tr1)';	Rk=forcesym(Rk);
    else
        Rk_inv=eye(ny)/Rk; Rk_inv=forcesym(Rk_inv);
    end
       
	%truncation 2
    if np>nm & strcmpi(doTruncation,'yes')
        [V2,lambda_orig2] = eig(J'*Rk_inv*J);
        [lambda2,I2] = sort(diag(lambda_orig2),1,'descend');
        lambdaJRkJ(:,k) = lambda2;
        V2 = V2(:,I2);
        tr2 = nm;
        JRkJ_inv= V2(:,1:tr2)*inv(diag(lambda2(1:tr2)))*V2(:,1:tr2)'; JRkJ_inv=forcesym(JRkJ_inv);
    else
        JRkJ_inv=eye(np)/(J'*Rk_inv*J); JRkJ_inv=forcesym(JRkJ_inv);   
    end
    
    Mk = JRkJ_inv*J'*Rk_inv;
%     p_filt(:,k)=Mk*(y(:,k)-G*x_pred(:,k));
    Ppkk=JRkJ_inv;  Ppkk=forcesym(Ppkk);
    
    %measurement update
    Lk=Pkk_*G'*Rk_inv;
%     x_filt(:,k)=x_pred(:,k)+Lk*(y(:,k)-G*x_pred(:,k)-J*p_filt(:,k));
    
    %truncation 3
	if ny>nm & strcmpi(doTruncation,'yes')  
        [V3,lambda_orig3] = eig(forcesym(J*Ppkk*J'));
        [lambda3,I3] = sort(diag(lambda_orig3),1,'descend');
        lambdaJPpkkJ(:,k) = lambda3;
        V3 = V3(:,I3);
        tr3 = min(nm,np);
        JPpkkJ = V3(:,1:tr3)*diag(lambda3(1:tr3))*V3(:,1:tr3)';  JPpkkJ=forcesym(JPpkkJ);
    else
        JPpkkJ=(J*Ppkk*J'); JPpkkJ=forcesym(JPpkkJ);   
    end
        
    Pkk = Pkk_ - Lk*(Rk-JPpkkJ)*Lk'; Pkk=forcesym(Pkk);
    Pxpkk=-Lk*J*Ppkk;
    Ppxkk=Pxpkk';

    %time update
%     x_pred(:,k+1)=A*x_filt(:,k)+B*p_filt(:,k);
    Nk=A*Lk*(eye(ny)-J*Mk)+B*Mk;
    Pkk_ = [A B]*[ Pkk Pxpkk;Ppxkk Ppkk]*[A';B']+Q-Nk*S'-S*Nk'; Pkk_=forcesym(Pkk_);

	% Trace
	Pkk_prev=Pkk_current; 	Pkk_current=Pkk;
    r_traceP(k)=trace(abs(Pkk_current-Pkk_prev))./trace(Pkk_current);

	Ppkk_prev=Ppkk_current; 	Ppkk_current=Ppkk;
    r_tracePp(k)=trace(abs(Ppkk_current-Ppkk_prev))./trace(Ppkk_current);
    
    if dispconv & mod(k,100)==0
        [r_traceP(k),r_tracePp(k)]
    end
    
    if k>minsteps & abs(r_traceP(k)) < convtol & abs(r_tracePp(k)) < convtol
		convreached=true;
		disp(['Trace convergence reached, k=' num2str(k)]);
		M_ss=Mk;
		L_ss=Lk;
		P_ss=Pkk;
		Pp_ss=Ppkk;
    elseif k>=maxsteps
    
% 		if strcmpi(doPlot,'yes')
		figure(); hold on; grid on; plot(r_traceP); set(gca,'YScale','log'); title('Ratio Px');
		figure(); hold on; grid on; plot(r_tracePp); set(gca,'YScale','log'); title('Ratio Pp');
% 		end
		error(['Divergence reached, k=' num2str(k)]);
    end
    
end

if strcmpi(doPlot,'yes')
figure(); hold on; grid on; plot(r_tracePp); set(gca,'YScale','log');
figure(); hold on; grid on; plot(r_traceP); set(gca,'YScale','log');
end

%% Filter estimates

for k=1:nt
    
    p_filt(:,k)=M_ss*(y(:,k)-G*x_pred(:,k));
    x_filt(:,k)=x_pred(:,k)+L_ss*(y(:,k)-G*x_pred(:,k)-J*p_filt(:,k));
    x_pred(:,k+1)=A*x_filt(:,k)+B*p_filt(:,k);
    
end

telapsed=toc(tstart);
disp(['JIS calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed*10^3./nt) ' seconds per 1k steps']);


if strcmpi(doScale,'yes')
  	x_filt=T1*x_filt;
  	P_ss=T1*P_ss*T1.';
	
	% Set these to empty for safety, not yet checked how these would be affected
  	M_ss=[];
  	L_ss=[];
end


%%













return
close all

figure(); hold on; grid on;
plot(lambdaRk'); xlim([0 k]);
set(gca,'YScale','log');


figure(); hold on; grid on;
plot(lambdaJRkJ');
set(gca,'YScale','log');

figure(); hold on; grid on;
plot(lambdaJPpkkJ');
set(gca,'YScale','log');

figure(); hold on; grid on;
plot(abs(lambdaJPpkkJ(1,:)./lambdaJPpkkJ(end,:)));
set(gca,'YScale','log');

%%
close all

figure(); hold on; grid on;
plot(squeeze(Rk_time(1,1,:)))

figure(); hold on; grid on;
plot(squeeze(Pkk_time(1,1,:)))

figure(); hold on; grid on;
plot(squeeze(Pkk__time(1,1,:)))

figure(); hold on; grid on;
plot(squeeze(Ppkk_time(1,1,:)))

set(gca,'YScale','log');

%%
close all

figure(); hold on; grid on;
plot(c1)
set(gca,'YScale','log');
title('Cony Rk')

figure(); hold on; grid on;
plot(c2)
set(gca,'YScale','log');
title('Cony J^T Rk J')

figure(); hold on; grid on;
plot(c3)
set(gca,'YScale','log');



