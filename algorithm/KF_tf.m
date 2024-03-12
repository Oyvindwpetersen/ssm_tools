function [H0y H1y H0p H1p]=KF_tf(A,B,G,J,Q,R,S,dt,omega_axis)
%% Transfer function for steady state operation of Kalman filter
%
% Model
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

[~,~,P_k_k,P_k_kmin]=KF(A,B,G,J,Q,R,S,y,p_dummy,[],[],'steadystate',true);

% [x_k_N,P_k_N]=RTSSmoother(Fad,x_k_k,x_k_kmin,P_k_k,P_k_kmin,'steadystate','yes');

% if any(any(S))
%     error('Implementation not supporting S~=0, due to modification necessary in RTS smoother');
% end

%% Calculate TF

P_k_kmin_ss=P_k_kmin;

Omega_k_ss=(G*P_k_kmin_ss*G.'+R); Omega_k_ss=forcesym(Omega_k_ss);
Omega_k_ss_inv=eye(size(Omega_k_ss))/Omega_k_ss; Omega_k_ss_inv=forcesym(Omega_k_ss_inv);
K_k_ss=(A*P_k_kmin_ss*G.'+S)*Omega_k_ss_inv;

M_k_ss=P_k_kmin_ss*G.'*Omega_k_ss_inv;

H_y=zeros(nx*2,ny,length(omega_axis));
H_p=zeros(nx*2,np,length(omega_axis));

matrix_y=[M_k_ss ; K_k_ss ];

for k=1:length(omega_axis)

    z=exp(1i*omega_axis(k)*dt);
    
    matrix_x=[
        eye(nx) , -eye(nx)+M_k_ss*G  ;
        zeros(nx) , z*eye(nx)-(A-K_k_ss*G)
        ];

    H_y(:,:,k)=matrix_x\matrix_y;
    
    if exist_input
        matrix_p=[-M_k_ss*J ; B-K_k_ss*J];
        H_p(:,:,k)=matrix_x\matrix_p;
    end

end

%% Split

H0y=H_y(1:nx,:,:);
H1y=H_y((nx+1):(nx*2),:,:);

if exist_input
    H0p=H_p(1:nx,:,:);
    H1p=H_p((nx+1):(nx*2),:,:);
else
    H0p=[];
    H1p=[];
end
