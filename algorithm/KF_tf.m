function [H_x0y_kf H_x1y_kf H_x0p_kf H_x1p_kf]=KF_tf(A,B,G,J,Q,R,S,dt,omega_axis)
%% Transfer function for steady state operation of Kalman filter
%
% Model:
% x(k+1)=A*x(k)+B*p(k)+w(k)
% y(k)=G*x(k)+J*p(k)+v(k)
%
% Inputs:
% A: state matrix
% B: input matrix
% G: output matrix
% J: direct transmission matrix
% Q: state noise covariance
% R: output noise covariance
% S: mixed noise covariance
% dt: time step
% omega_axis: frequency axis
%
% Outputs:
% H0y: transfer from y to filter estimate
% H1y: transfer from y to prediction estimate
% H0p: transfer from p to filter estimate
% H1p: transfer from p to prediction estimate
%
%
% | x0(w) | = | H0y(w) | y(w) + | H0p(w) | p(w)
% | x1(w) |   | H1y(w) |        | H1p(w) | 
%

%% Parameters

nx=size(A,1);
ny=size(G,1);
np=size(B,2);

if isempty(B) & isempty(J)
    exist_input=false;
    p_dummy=[];
else
    exist_input=true;
    p_dummy=zeros(size(B,2),1e3);
end

y=zeros(size(G,1),1e3);

%% Run KF

[~,~,P_k_k,P_k_kmin,K_k_ss,M_k_ss]=KF(A,B,G,J,Q,R,S,y,p_dummy,[],[],'steadystate',true);

%% Calculate TF

H_y=zeros(nx*2,ny,length(omega_axis));
H_p=zeros(nx*2,np,length(omega_axis));

matrix_y=[M_k_ss ; K_k_ss ];

if exist_input
    matrix_p=[-M_k_ss*J ; B-K_k_ss*J];
end

for k=1:length(omega_axis)

    z=exp(1i*omega_axis(k)*dt);
    
    matrix_x=[
        eye(nx) , -eye(nx)+M_k_ss*G  ;
        zeros(nx) , z*eye(nx)-(A-K_k_ss*G)
        ];

    H_y(:,:,k)=matrix_x\matrix_y;
    
    if exist_input
        H_p(:,:,k)=matrix_x\matrix_p;
    end

end

%% Split

H_x0y_kf=H_y(1:nx,:,:);
H_x1y_kf=H_y((nx+1):(nx*2),:,:);

if exist_input
    H_x0p_kf=H_p(1:nx,:,:);
    H_x1p_kf=H_p((nx+1):(nx*2),:,:);
else
    H_x0p_kf=[];
    H_x1p_kf=[];
end
