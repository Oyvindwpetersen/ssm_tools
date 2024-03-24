function [F,L,H,Ht,S_par,R_par]=ssmod_maternglobal(lambda_vec,sigma_w_vec,p_vec,omega_axis,tau_axis,calcfast)

%% Block diagonal Matern state space model
%
% Inputs:
% lambda_vec: vector with hyperparameters, inverse length scale
% sigma_vec: vector with hyperparameters, standard deviations
% p_vec: vector with order
% omega_axis: frequency axis in [rad/s], can be omitted
% tau_axis: time axis in [s], can be omitted
% calcfast: false/true for premade equations or general calculation through state-space model
%
% Outputs:
% F: state matrix (cont)
% L: input matrix (cont)
% H: output matrix (cont)
%
%%

if nargin<6
    calcfast=false;
end

if nargin<5
    tau_axis=[];
end

if nargin<4
    omega_axis=[];
end

if any(sigma_w_vec<0)
    error('Only positive sigma allowed');
end

if any(lambda_vec<0)
    error('Only positive lambda allowed');
end

%%

n_comp=length(sigma_w_vec);

F_cell=cell(1,n_comp);
L_cell=cell(1,n_comp);
H_cell=cell(1,n_comp);

if calcfast
	S_par=zeros(n_comp,n_comp,length(omega_axis));
end

for k=1:n_comp
    
	lambda=lambda_vec(k);
	sigma_w=sigma_w_vec(k);

	[F_cell{k},L_cell{k},H_cell{k}]=ssmod_matern(lambda,p_vec(k),[]);

	% Two-sided PSD
	if calcfast
		S_par(k,k,:)=sum( 1/(2*pi)*(sigma_w.^2) .* ( (lambda.^2)+omega_axis.^2 ).^(-p_vec(k)-1) , 1);
		Ht=[];
	end
end

% CF
if calcfast
    R_par=matern_cf(tau_axis,lambda_vec,sigma_w_vec,p_vec,'3d');
end

F=blkdiag(F_cell{:});
L=blkdiag(L_cell{:});
H=blkdiag(H_cell{:});

%% Slow way

if ~calcfast

	% TF for state-space model
	Ht=ssmod_tf(F,L,H,0,omega_axis,[],'type','io');

	% PSD for white noise
	Sigma_w_squared=repmat(diag(sigma_w_vec).^2,1,1,length(omega_axis));
	S_w=1/(2*pi)*Sigma_w_squared;

	% Two-sided PSD
	S_par=mtimes3(Ht,S_w,Ht,'nnh');

	% CF
	Sigma_w_squared=diag(sigma_w_vec).^2;
	Q=L*Sigma_w_squared*L.';
	[Pinf,~,~] = icare(F.',[],Q,[],[],[],[]);
	phiFunc=@(tau) expmnorm(F*tau);
	R_par=zeros(size(H,1),size(H,1),length(tau_axis));
	for k=1:length(tau_axis)
		if tau_axis(k)<0;
			error('Only positive tau allowed');
		end
		R_par(:,:,k)=H*Pinf*phiFunc(tau_axis(k)).'*H.';
	end

end
