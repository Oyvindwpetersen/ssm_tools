function [x_filt p_filt x_pred]=DKF(A,B,G,J,y,x0,p0,Q,R,Qp,P0,Pp0,varargin)

%% Dual Kalman Filter
%
% Inputs:
% A: state matrix
% B: input matrix
% G: output matrix
% J: direct transmission matrix
% y: output vector
% x0: initial state
% p0: initial input
% R: output noise covariance
% Q: state noise covariance
% P0: initial state error covariance
% Pp0: initial input error covariance
%
% Outputs:
% x_filt: filter state estimate
% p_filt: filter input estimate
% x_pred: prediction state estimate

%% Inputs
p=inputParser;
addParameter(p,'steadystate','yes',@ischar)
parse(p,varargin{:});
steadystate = p.Results.steadystate;

%% Initiate

ns=size(A,1);
ny=size(G,1);
np=size(B,2);
nt=length(y);

% zero matrices
x_pred=zeros(ns,nt);
x_filt=zeros(ns,nt);
p_pred=zeros(np,nt);
p_filt=zeros(np,nt);

% assign initial values
x_filtk_=[x0];
p_filtk_=[p0];
Pk_=P0;
Ppk_=Pp0;

%%

if strcmpi(steadystate,'no')

tstart=tic;
for k=1:nt

	%input
	p_pred(:,k)=p_filtk_;
	Pp_k=Ppk_+Qp;  Pp_k=forcesym(Pp_k);
	
	%kalman for input
    
	Gpk1=Pp_k*J' / (J*Pp_k*J' +R);
	
	p_filt(:,k)=p_pred(:,k)+Gpk1*( y(:,k) - G*x_filtk_ - J*p_pred(:,k) );
	Ppk=Pp_k-Gpk1*J*Pp_k; Ppk=forcesym(Ppk);
	
	%state
	x_pred(:,k)=A*x_filtk_+B*p_filt(:,k);
	P_k=A*Pk_*A'+Q; P_k=forcesym(P_k);
	
	%kalman for state
	Gxk2=P_k*G' / (G*P_k*G' + R);	
	x_filt(:,k)=x_pred(:,k)+Gxk2*( y(:,k) - G*x_pred(:,k) - J*p_filt(:,k) );
    Pk=P_k-Gxk2*G*P_k; Pk=forcesym(Pk);
   	
	%vectors and matrices for next step: k -> k-1
	p_filtk_=p_filt(:,k);
	x_filtk_=x_filt(:,k);
    Ppk_=Ppk;
    Pk_=Pk;

end
telapsed=toc(tstart);

disp(['DKF calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed*10^3./nt) ' seconds per 1k steps']);

end

%%

if strcmpi(steadystate,'yes')
    
nt_ss=1000;
nt_check=100;

trace_Pk=zeros(1,nt_ss);
trace_Ppk=zeros(1,nt_ss);

diag_Pk=zeros(ns,nt_check);
diag_Ppk=zeros(np,nt_check);

tstart=tic;

for k=1:nt_ss

	%input
	Pp_k=Ppk_+Qp;  Pp_k=forcesym(Pp_k);
	
	%kalman for input
	Gpk1=Pp_k*J' / (J*Pp_k*J' +R);
	
	Ppk=Pp_k-Gpk1*J*Pp_k; Ppk=forcesym(Ppk);
	
	%state
	P_k=A*Pk_*A'+Q; P_k=forcesym(P_k);
	
	%kalman for state
	Gxk2=P_k*G' / (G*P_k*G' + R);	
    Pk=P_k-Gxk2*G*P_k; Pk=forcesym(Pk);
   	
	%vectors and matrices for next step: k -> k-1
    Ppk_=Ppk;
    Pk_=Pk;
    
    % trace
    trace_Pk(k)=trace(Pk);
    trace_Ppk(k)=trace(Ppk);
    
    diag_Pk=[ diag_Pk(:,2:nt_check) diag(Pk)];
    diag_Ppk=[ diag_Ppk(:,2:nt_check) diag(Ppk)];    
    
end

trace_r1=max( abs(trace_Pk(end)./trace_Pk(end-nt_check) -1));
trace_r2=max( abs(trace_Ppk(end)./trace_Ppk(end-nt_check) -1));
if trace_r1>1e-3 | trace_r2>1e-3
    warning(['Trace ratio state: ' num2str(trace_r1,'%.5e') ', trace ratio force: ' num2str(trace_r2,'%.5e')]);
end

diag_r1= max(max(abs( diag_Pk ./ repmat(diag_Pk(:,end),1,nt_check)-1),[],2));
diag_r2= max(max(abs( diag_Ppk ./ repmat(diag_Ppk(:,end),1,nt_check)-1),[],2));
if diag_r1>1e-3 | diag_r2>1e-3
    warning(['Ratio state: ' num2str(diag_r1,'%.5e') ', ratio force: ' num2str(diag_r2,'%.5e')]);
end
    
for k=1:nt

	%input
	p_pred(:,k)=p_filtk_;
% 	Pp_k=Ppk_+Qp;  Pp_k=forcesym(Pp_k);
	
	%kalman for input
% 	Gpk=Pp_k*J' / (J*Pp_k*J' +R);
	
	p_filt(:,k)=p_pred(:,k)+Gpk1*( y(:,k) - G*x_filtk_ - J*p_pred(:,k) );
% 	Ppk=Pp_k-Gpk*J*Pp_k; Ppk=forcesym(Ppk);
	
	%state
	x_pred(:,k)=A*x_filtk_+B*p_filt(:,k);
% 	P_k=A*Pk_*A'+Q; P_k=forcesym(P_k);
	
	%kalman for state
% 	Gxk=P_k*G' / (G*P_k*G' + R);	
	x_filt(:,k)=x_pred(:,k)+Gxk2*( y(:,k) - G*x_pred(:,k) - J*p_filt(:,k) );
%     Pk=P_k-Gxk*G*P_k; Pk=forcesym(Pk);
   	
	%vectors and matrices for next step: k -> k-1
	p_filtk_=p_filt(:,k);
	x_filtk_=x_filt(:,k);
%     Ppk_=Ppk;
%     Pk_=Pk;

end
    
telapsed=toc(tstart);
disp(['DKF calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed*10^3./nt) ' seconds per 1k steps']);
    
end

