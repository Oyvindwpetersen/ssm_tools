function [A_mat,B,R_target,tau_target]=AR_fit_old(w_axis,S_target,j_lag,l_lag,dt)

%% Fit AR-model with sparse lags
% From Gallego-Castillo 2021, A tutorial on reproducing a predefined autocovariance function
% through AR models: application to stationary homogeneous isotropic turbulence

% z[n] = sum_r_1toP( A[j(r)]*z[n-j(r)] ) + B*epsilon[n]
% Cov(epsilon)=E[epsilon*epsilon^T]=I

% Inputs:
% w_axis: frequency axis in [rad/s]
% S: spectral density
% j_lag: lags in AR-model, can be regular [1,2,3] or skipping [1,3,9]
% l_lag: lags in AR-model (target)
% dt: time discretization

% Outputs:
% S_AR: spectral density of AR-model
% H: transfer function of AR-model
% A_mat: 3d-form of A-matrices
% B: input matrix
% R: 
% tau: 

%%

j_lag=j_lag(:);
l_lag=l_lag(:);

% lag_l_j=repmat(l_lag.',length(l_lag),1)-repmat(j_lag,1,length(j_lag));

for ind1=1:length(j_lag)
for ind2=1:length(l_lag)
    lag_l_j(ind1,ind2)=l_lag(ind2)-j_lag(ind1);
end
end

lag_calc=unique( [ 0 ; lag_l_j(:) ; l_lag(:) ; j_lag(:) ]);
lag_calc=[min(lag_calc):max(lag_calc)];
tau_target=dt*lag_calc;

R_target=zeros(size(S_target,1),size(S_target,1),length(lag_calc));

for k=length(lag_calc):-1:1
    
    if lag_calc(k)<0
        
        ind=find(abs(lag_calc(k))==lag_calc);
        R_target(:,:,k)=R_target(:,:,ind).';
        
    else
    
        cos_vec(1,1,:)=cos(w_axis*dt*lag_calc(k));
        cos_mat=repmat(cos_vec,size(S_target,1),size(S_target,1),1);

        integrand=S_target.*cos_mat;

    %     R_target(:,:,k)=trapz(w_axis,integrand,3);
        R_target(:,:,k)=simps(w_axis,integrand,3);

    end
    
end


%% Fit AR model

Mat_l_j=cell(length(j_lag),length(j_lag));
for ind1=1:length(j_lag)
    for ind2=1:length(l_lag)
        j=lag_l_j(ind1,ind2);
        ind3=find(lag_calc==j);
        Mat_l_j{ind1,ind2}=R_target(:,:,ind3);
    end
end
Mat_l_j=cell2mat(Mat_l_j);


Mat_l=cell(1,length(l_lag));
for ind1=1:length(l_lag)
    j=l_lag(ind1);
    ind3=find(lag_calc==j);
    Mat_l{1,ind1}=R_target(:,:,ind3);
end
Mat_l=cell2mat(Mat_l);


Mat_j=cell(1,length(j_lag));
for ind1=1:length(j_lag)
    j=j_lag(ind1);
    ind3=find(lag_calc==j);
    Mat_j{1,ind1}=R_target(:,:,ind3);
end
Mat_j=cell2mat(Mat_j);


% Find A_block=[A1,A2,A3,...]
A_block=Mat_l/Mat_l_j;

A_mat=reshape(A_block,size(A_block,1),size(A_block,1),length(j_lag));

ind3=find(lag_calc==0);
Mat0=R_target(:,:,ind3);

BBt=Mat0-A_block*Mat_j.';

B=chol(BBt);

% R=R_target;

% [S_AR,H_AR]=AR_PSD(w_axis,A_mat,B,j_lag,dt);

