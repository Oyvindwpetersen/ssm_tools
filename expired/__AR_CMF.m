function [S_AR,H_AR,A_mat,B,R,tau]=AR_CMF(w_axis,S,j_lag,dt)

%% Fit AR-model with sparse lags
% From Møller 2019, Turbulent wind field representation and conditional mean field simulation

% y[n] = sum_r_1toP( A_r[j(r)]*y[n-j(r)] ) + B*eta[n]

% Cov(eta)=E[eta*eta^T]=I

% Inputs:
% w_axis: frequency axis in [rad/s]
% S: spectral density
% j_lag: lags in AR-model, can be regular [1,2,3] or skipping [1,3,9]
% dt: time discretization

% Outputs:
% S_AR: spectral density of AR-model
% H: transfer function of AR-model
% A_mat: 3d-form of A-matrices
% B: input matrix
% R: 
% tau: 

%%

% Cww
j_lag_comb=repmat(-j_lag.',1,length(j_lag))+repmat(j_lag,length(j_lag),1);

j_lag_comb_unique=unique( [ j_lag_comb(:) ; ; max(j_lag) ; -max(j_lag) ; j_lag(:) ]);

R_lag_comb_unique=zeros(size(S,1),size(S,1),length(j_lag_comb_unique));

for k=length(j_lag_comb_unique):-1:1
    
    if j_lag_comb_unique(k)<0
        
        ind=find(abs(j_lag_comb_unique(k))==j_lag_comb_unique);
        R_lag_comb_unique(:,:,k)=R_lag_comb_unique(:,:,ind).';
        
    else
    
        cos_vec(1,1,:)=cos(w_axis*dt*j_lag_comb_unique(k));
        cos_mat=repmat(cos_vec,size(S,1),size(S,1),1);

        integrand=S.*cos_mat;

    %     R_lag_comb_unique(:,:,k)=trapz(w_axis,integrand,3);
        R_lag_comb_unique(:,:,k)=simps(w_axis,integrand,3);

    end
    
end


%% Solve AR model

Cww_cell=cell(length(j_lag),length(j_lag));
for ind1=1:length(j_lag)
for ind2=1:length(j_lag)
    
    j=j_lag_comb(ind1,ind2);
    ind3=find(j_lag_comb_unique==j);
    Cww_cell{ind1,ind2}=R_lag_comb_unique(:,:,ind3);
    
end
end

Cww=cell2mat(Cww_cell);

Cuw_cell=cell(1,length(j_lag));
for ind1=1:length(j_lag)
    j=j_lag(ind1);
    ind3=find(j_lag_comb_unique==j);
    Cuw_cell{1,ind1}=R_lag_comb_unique(:,:,ind3);
    
end
Cuw=cell2mat(Cuw_cell);

ind3=find(j_lag_comb_unique==0);

Cuu=R_lag_comb_unique(:,:,ind3);

A_block=Cuw/Cww;
% A_block=(Cww\(Cuw.')).';
A_mat=reshape(A_block,size(A_block,1),size(A_block,1),length(j_lag));

BBt=Cuu-A_block*Cww*A_block.';

B=chol(BBt);

[S_AR,H_AR]=AR_PSD(w_axis,A_mat,B,j_lag,dt);

