function [negLL,negLL_1,negLL_2,negLL_1h,negLL_2h,Nred]=loglik_evidence(H,x_k_kmin,y,Omega,varargin)

%%

p=inputParser;

addParameter(p,'cut',[0 0],@isnumeric)
parse(p,varargin{:});

cut=p.Results.cut;


%% Log likelihood of model evidence in state-space function

% Inputs:
% H: Output matrix
% x_k_kmin: prediction estimate from Kalman filter
% y: output data
% Omega: covariance matrix of innovation, Omega=(H*Px*H^T+R)

% p(y|theta)

% negLL= sum_k_1_N (  log(det(Omega_k)) + e_k^T*inv(Omega)*e_k )

%%
if cut(1)>0 | cut(2)>0;
    
    ind_start=[1+cut(1)];
    ind_end=size(y,2)-cut(2);
    
    range_keep=[ind_start:ind_end];
    y=y(:,range_keep);
    x_k_kmin=x_k_kmin(:,range_keep);
    
    Nred=size(y,2);
    
else
    Nred=size(y,2);
end

%% Model complexity

N=size(y,2);

negLL_1=N*log(det(Omega));

negLL_1h=negLL_1/N;

%% Data fit

e=y-H*x_k_kmin;

Omega=forcesym(Omega);
[V,D]=eig(Omega);

Omega_half=V*D.^0.5;
Omega_half_inv=eye(size(Omega_half))/Omega_half;
Temp2=Omega_half_inv*e;
negLL_2_temp=sum(Temp2.^2,1);
negLL_2=sum(negLL_2_temp);

negLL_2h=negLL_2_temp;

% negLL_2_check=0;
% for k=1:length(e)
%     negLL_2_check=negLL_2_check+e(:,k).'/Omega*e(:,k);
% end

%% Sum

negLL=negLL_1+negLL_2;

%% Test code

% no=30
% nt=12e3*20
% 
% Omega=randn(no); Omega=Omega*Omega.'+eye(no);
% e=randn(no,nt);
% y=nan*e;


% 
% tic
% Temp1=0;
% for k=1:length(e)
%     Temp1=Temp1+e(:,k).'/S*e(:,k);
% end
% toc
% 
% tic
% [V,D]=eig(S);
% 
% S_half=V*D.^0.5;
% S_half_inv=eye(size(S_half))/S_half;
% 
% Temp2=S_half_inv*e;
% 
% Temp3=sum(sum(Temp2.^2,1),2);
% toc
% 
% Temp1
% Temp3