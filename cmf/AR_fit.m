function [A_mat,B]=AR_fit(tau_target,R_target,j_lag,l_lag,dt)

%% Fit AR-model with sparse lags
% From Gallego-Castillo 2021, A tutorial on reproducing a predefined autocovariance function
% through AR models: application to stationary homogeneous isotropic turbulence
%
% z[n] = sum_r_1toP( A[j(r)]*z[n-j(r)] ) + B*epsilon[n]
% Cov(epsilon)=E[epsilon*epsilon^T]=I
%
% Inputs:
% tau_target: frequency axis in [rad/s]
% R_target: spectral density
% j_lag: lags in AR-model, can be regular [1,2,3] or skipping [1,3,9]
% l_lag: lags in AR-model (target in fit)
% dt: time discretization
%
% Outputs:
% A_mat: 3d-form of A-matrices
% B: input matrix
%
%%

j_lag=j_lag(:);
l_lag=l_lag(:);

for ind1=1:length(j_lag)
for ind2=1:length(l_lag)
    lag_l_j(ind1,ind2)=l_lag(ind2)-j_lag(ind1);
end
end

lag=tau_target/dt;

lag_rounded=round(lag);

if any(abs(lag_rounded-lag)>1e-6)
    error('Lags must coincide with integer steps of dt')
else
    lag=lag_rounded;
end

%% Fit AR model

Mat_l_j=cell(length(j_lag),length(j_lag));
for ind1=1:length(j_lag)
    for ind2=1:length(l_lag)
        j=lag_l_j(ind1,ind2);
        
        if j<0
            j=abs(j);
            ind3=find(lag==j);
            Mat_l_j{ind1,ind2}=R_target(:,:,ind3).';
        else
            ind3=find(lag==j);
            Mat_l_j{ind1,ind2}=R_target(:,:,ind3);
        end
        
        if isempty(ind3)
            j
            error('Missing this time lag');
            
        end
    end
end
Mat_l_j=cell2mat(Mat_l_j);


Mat_l=cell(1,length(l_lag));
for ind1=1:length(l_lag)
    j=l_lag(ind1);
    ind3=find(lag==j);
    Mat_l{1,ind1}=R_target(:,:,ind3);
end
Mat_l=cell2mat(Mat_l);


Mat_j=cell(1,length(j_lag));
for ind1=1:length(j_lag)
    j=j_lag(ind1);
    ind3=find(lag==j);
    Mat_j{1,ind1}=R_target(:,:,ind3);
end
Mat_j=cell2mat(Mat_j);


% Find A_block=[A1,A2,A3,...]
A_block=Mat_l/Mat_l_j;

A_mat=reshape(A_block,size(A_block,1),size(A_block,1),length(j_lag));

ind3=find(lag==0);
Mat0=R_target(:,:,ind3);

BBt=Mat0-A_block*Mat_j.';

B=chol(BBt);

% R=R_target;

% [S_AR,H_AR]=AR_PSD(w_axis,A_mat,B,j_lag,dt);

