function [x_filt p_filt P_ss Pp_ss M_ss L_ss] = JIS_trunc_ss2(A,B,G,J,y,x0,R,Q,S,P01,varargin)

%% Implementation from Niu (2011)
% Taking into account S


%% Input

p=inputParser;
addParameter(p,'trunc','no',@ischar)
addParameter(p,'plot','no',@ischar)
addParameter(p,'scale','no',@ischar)

parse(p,varargin{:});
doTruncation = p.Results.trunc;
doPlot = p.Results.plot;
doScale = p.Results.scale;

%%


%%

R_inv=eye(size(R))/R;

A_star=A-S*R_inv*G;
B_star=B-S*R_inv*J;

Q_star=Q-S*R_inv*S.';
S_star=zeros(size(S));
R_star=R;
clear A B Q %R S

%%

ns=size(A_star,1);
nt=length(y);
nd=size(G,1);
np=size(B_star,2);
nm=ns/2;
minTimeSteps=100;
maxTimeSteps=100e3;
convTol=1e-6;

%zero matrices
x_pred=zeros(ns,nt);
x_filt=zeros(ns,nt);
p_filt=zeros(np,nt);
% innov=zeros(nd,nt);

lambdaRk=zeros(nd,maxTimeSteps);
lambdaJRkJ=zeros(np,maxTimeSteps);
lambdaJPpkkJ=zeros(nd,maxTimeSteps);

%assign initial values
x_pred(:,1)=[x0];
Pkk_=P01;


% [y G J R S]=scaleGJ(y,G,J,R,S,diag(R));
% [A,B,G,J,y,Q,R,S,T1,T2]=scaleSysUnitCov(A,B,G,J,y,Q,R,S);
% [~,~,G,J,y,~,R,S,~,~]=scaleSysUnitCov([],[],G,J,y,[],R,S);

% if strcmpi(doScale,'yes')
   % [A,B,G,J,y,Q,R,S,T1,T2]=scaleSysUnitCov(A,B,G,J,y,Q,R,S); 
% end

%% Steady state

% trace_Pkk=zeros(1,nt_ss);
% trace_Ppkk=zeros(1,nt_ss);

clear r_traceP r_tracePp

Pkk_prev=zeros(ns); Pkk_current=zeros(ns);
Ppkk_prev=zeros(np); Ppkk_current=zeros(np);

tstart=tic;
convReached=false; k=0;
while convReached==false    
    k=k+1;
    
    %input estimation
    Rk =G*Pkk_*G'+R_star;  Rk=forcesym(Rk); 
%     nd
%     nm
    %truncation 1
    if 0 %nd>nm% & strcmpi(doTruncation,'yes');  %%%%%%%%%%%%%%%%%%%%%%%%%%%% - disenabled on BSB for reasons of ?
        [V1,lambda_orig1] = eig(Rk);
        [lambda1,I1] = sort(diag(lambda_orig1),1,'descend');
        lambdaRk(:,k) = lambda1;
        V1 = V1(:,I1);
        tr1 = nm;
        Rk_inv = V1(:,1:tr1)*inv(diag(lambda1(1:tr1)))*V1(:,1:tr1)'; Rk_inv=forcesym(Rk_inv);
        Rk = V1(:,1:tr1)*diag(lambda1(1:tr1))*V1(:,1:tr1)';	Rk=forcesym(Rk);
    else
        Rk_inv=eye(nd)/Rk; Rk_inv=forcesym(Rk_inv);
    end
       
	%truncation 2
    if np>nm & strcmpi(doTruncation,'yes'); 
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
	if nd>nm & strcmpi(doTruncation,'yes');  
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
    % Nk=A*Lk*(eye(nd)-J*Mk)+B*Mk;
    Pkk_ = [A_star B_star]*[ Pkk Pxpkk;Ppxkk Ppkk]*[A_star';B_star']+Q_star;%-Nk*S'-S*Nk'; Pkk_=forcesym(Pkk_);

	% trace
        
	Pkk_prev=Pkk_current; 	Pkk_current=Pkk;
    r_traceP(k)=trace(abs(Pkk_current-Pkk_prev))./trace(Pkk_current);

	Ppkk_prev=Ppkk_current; 	Ppkk_current=Ppkk;
    r_tracePp(k)=trace(abs(Ppkk_current-Ppkk_prev))./trace(Ppkk_current);
    
    if k>minTimeSteps & abs(r_traceP(k)) < convTol & abs(r_tracePp(k)) < convTol
    convReached=true;
    disp(['Trace convergence reached, k=' num2str(k)]);
    M_ss=Mk;
    L_ss=Lk;
    P_ss=Pkk;
    Pp_ss=Ppkk;
    elseif k>maxTimeSteps
    
    if strcmpi(doPlot,'yes')
    figure(); hold on; grid on; plot(r_tracePp); set(gca,'YScale','log');
    figure(); hold on; grid on; plot(r_traceP); set(gca,'YScale','log');
    end
    error(['Divergence reached, k=' num2str(k)]);
    end
    
%     if mod(k,10)==0
%     [r_traceP(k) r_tracePp(k)]
%     end

end

if strcmpi(doPlot,'yes')
figure(); hold on; grid on; plot(r_tracePp); set(gca,'YScale','log');
figure(); hold on; grid on; plot(r_traceP); set(gca,'YScale','log');
end

%% Filter estimates

for k=1:nt
    
    p_filt(:,k)=M_ss*(y(:,k)-G*x_pred(:,k));
    x_filt(:,k)=x_pred(:,k)+L_ss*(y(:,k)-G*x_pred(:,k)-J*p_filt(:,k));
    x_pred(:,k+1)=A_star*x_filt(:,k)+B_star*p_filt(:,k)+S*R_inv*y(:,k);
    
end

telapsed=toc(tstart);
disp(['JIS calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed*10^3./nt) ' seconds per 1k steps']);




