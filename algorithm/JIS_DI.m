function [x_filt p_filt Px Pp]=JIS_DI(A,B,B_prime,C,D,D_prime,y,u,x0,R,Q,S,P01,varargin)

%% Joint input and state estimation for system with additional known inputs:
% x(k+1)=A*x(k)+B*p(k)+B_prime*u(k)+w(k);
% y=C*x(k)+D*p(k)+D_prime*u(k)+v(k);

% From Song (2018)

% Inputs:
% A: state matrix
% B: input matrix
% B_prime: input matrix for known input
% C: output matrix
% D: direct transmission matrix
% D_prime: direct transmission matrix for known input
% y: output vector
% u: known input
% x0: initial state estimate
% R: output noise covariance
% Q: state noise covariance
% S: mixed noise covariance
% P01: initial state error covariance

% Outputs:
% x_filt: filter state estimate
% p_filt: filter input estimate
% Px: filter state error covariance
% Pp: filter input error covariance

p=inputParser;
addParameter(p,'steadystate','no',@ischar)
parse(p,varargin{:});
steadystate = p.Results.steadystate;

%% Initiate

ns=size(A,1);
nt=length(y);
np=size(B,2);

% Zero matrices
x_k_k=zeros(ns,nt);
x_k_kmin=zeros(ns,nt);
p_k_k=zeros(np,nt);

Px_k_k=zeros(ns,ns,nt);
Pp_k_k=zeros(np,np,nt);

% Assign initial values
x_k_kmin(:,1)=[x0];
Px_k_kmin=P01;

%% Conventional one step recursive

if strcmpi(steadystate,'no');

tstart=tic;
for k=1:nt;
	
	% 
	if k>1
		Px_k_kmin=Px_kplus_k;
	else
		Px_k_kmin=P01;
	end
	
	%%% Time update
	% Eq B3 and B4 in Song
	% x_k_kmin(:,k)=A*x_k_k(:,k-1)+B_prime*u(:,k-1)+B*p_k_k(:,k-1);
	% Px_k_kmin=[A B]*[Px_kmin_kmin Pxp_kmin_kmin ; Ppx_kmin_kmin Pp_kmin_kmin]*[A.' ; B.']+Q

	%%% Input estimation
	% Eq B6 in Song
	Pe_k=C*Px_k_kmin*C.'+R;
    
	% Eq 13 in Song
	Temp=D.'/Pe_k*D;
	M_k=Temp\D.'/Pe_k;
	
	% Just before Eq 8 in Song
	% y_tilde is the innovation
	y_tilde(:,k)=y(:,k)-C*x_k_kmin(:,k)-D_prime*u(:,k);
	
	p_k_k(:,k)=M_k*y_tilde(:,k);

    Pp_k_k=M_k*Pe_k*M_k.';
    
	%%% Measurement update
	
	% B7 in Song
	Pxe_k_kmin=Px_k_kmin*C.';
	Pex_k_kmin=Pxe_k_kmin.';

	% Page 892 in Song
	K_k=Pxe_k_kmin/Pe_k;
	
	L_k=K_k*(eye(size(D,1))-D*M_k);
	
	x_k_k(:,k)=x_k_kmin(:,k)+K_k*(eye(size(D,1))-D*M_k)*(y(:,k)-C*x_k_kmin(:,k)-D_prime*u(:,k));
	
	Px_k_k=Px_k_kmin-K_k*(Pe_k-D*Pp_k_k*D.')*K_k.';
	
	Pxp_k_k=-Pxe_k_kmin*M_k.';
	Ppx_k_k=Pxp_k_k.';
    
	%%% Time update
	% Eq B3 and B4 in Song
	x_k_kmin(:,k+1)=A*x_k_k(:,k)+B_prime*u(:,k)+B*p_k_k(:,k);
	
	Px_kplus_k=[A B]*[Px_k_k Pxp_k_k ; Ppx_k_k Pp_k_k]*[A.' ; B.']+Q;
	
	% Save covariance matrices
	Px_k_k_all(:,:,k)=Px_k_k;
	Pp_k_k_all(:,:,k)=Pp_k_k;
	
end
telapsed=toc(tstart);
disp(['JIS calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed*10^3./nt) ' seconds per 1k steps']);

end

%% Output

x_filt=x_k_k;
p_filt=p_k_k;
Px=Px_k_k_all;
Pp=Pp_k_k_all;


