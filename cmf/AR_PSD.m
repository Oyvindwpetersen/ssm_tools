function [S,H]=AR_PSD(w_axis,A,B,j_lag,dt)

%% PSD of AR-model

% z[n] = sum_r_1toP( A[j(r)]*z[n-j(r)] ) + B*epsilon[n]
% Cov(epsilon)=E[epsilon*epsilon^T]=I

% Inputs:
% w_axis: frequency axis in [rad/s]
% A: block of square A-matrices: A=[A_1,A_2,A_3,...] or in 3d: A(:,:,1)=A_1, A(:,:,2)=A_2, A(:,:,3)=A_3
% B: input matrix
% j_lag: lags in AR-model, can be regular [1,2,3] or skipping [1,3,9]
% dt: time lag

% Outputs:
% S: spectral density of output
% H: transfer function of AR-model





% The standard deviation of epsilon (in disc time) is 1. 
% In general, Sigma_eps_disc^2=Sigma_eps_cont^2*dt
% The PSD of of epsilon is defined as 

% S_eta=Sigma_eps_cont^2/(2*pi)=



%%

A_mat=AR_coeffconv(A,'matrix');

H=zeros(size(A_mat,1),size(A_mat,2),length(w_axis));

for k=1:length(w_axis)
    
    A_exp=zeros(size(A_mat));
    for r=1:size(A_mat,3)
        A_exp(:,:,r)=-A_mat(:,:,r)*exp(-1i*w_axis(k)*j_lag(r)*dt);
    end
    
    H(:,:,k)= ( eye(size(A_mat,1)) / ( eye(size(A_mat,1)) + sum(A_exp,3)) ); %*exp(-w_axis(k)*r*dt);

end

S_eta=repmat(B*B.',1,1,length(w_axis))*1/(2*pi)*dt;

% Todo: find out why dt is needed here. Should not be needed in already disc models?

S=mtimes3(H,S_eta,H,'nnh');

%%



