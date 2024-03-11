function [H0y H1y H0p H1p]=KF_tf(A,B,G,J,Q,R,S,dt,omega_axis)

%%

%% Run KF

% x0=zeros(size(Fad,1),1);
% P_0_0=[];

y=zeros(size(G,1),1e3);
p_det=zeros(size(B,2),1e3);


if isempty(B) & isempty(J)
    exist_input=false;
else
    exist_input=true;
end   

% Obtain matrices

[x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KF(A,B,G,J,Q,R,S,y,p_det,[],[],'steadystate',true);

% [x_k_N,P_k_N]=RTSSmoother(Fad,x_k_k,x_k_kmin,P_k_k,P_k_kmin,'steadystate','yes');

% if any(any(S))
%     error('Implementation not supporting S~=0, due to modification necessary in RTS smoother');
% end

%% TF of filter and smoother

P_k_kmin_ss=P_k_kmin;
P_k_k_ss=P_k_k;

Omega_k_ss=(G*P_k_kmin_ss*G.'+R); Omega_k_ss=forcesym(Omega_k_ss);
Omega_k_ss_inv=eye(size(Omega_k_ss))/Omega_k_ss; Omega_k_ss_inv=forcesym(Omega_k_ss_inv);
K_k_ss=(A*P_k_kmin_ss*G.'+S)*Omega_k_ss_inv;

M_k_ss=P_k_kmin_ss*G.'*Omega_k_ss_inv;

nx=size(A,1);
ny=size(G,1);
np=size(B,2);

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

H0y=H_y(1:nx,:,:);
H1y=H_y((nx+1):(nx*2),:,:);

if exist_input
    H0p=H_p(1:nx,:,:);
    H1p=H_p((nx+1):(nx*2),:,:);
else
    H0p=[];
    H1p=[];
end


% for k=1:length(omega_axis)
%     z=exp(1i*w_axis(k)*mod.dt);
%     Mat=[
%         -eye(nx)+K_k_ss*mod.G eye(nx) zeros(nx) ;
%         z*eye(nx) -mod.A zeros(nx) ;
%         N_k_ss*z -eye(nx) eye(nx)-N_k_ss*z
%         ];
%
%     H_yx(:,:,k)=Mat\[K_k_ss ; zeros(nx,ny) ; zeros(nx,ny) ];

% end



%
% for k=1:length(w_axis)
%     z=exp(1i*w_axis(k)*dt);

%     Mat=[
%         -eye(nx)+K_k_ss*Had eye(nx) zeros(nx) ;
%         z*eye(nx) -Fad zeros(nx) ;
%         N_k_ss*z -eye(nx) eye(nx)-N_k_ss*z
%         ];

%     Mat_L=[
%         -eye(nx)+P_k_kmin_ss*Had.'*Omega_k_ss_inv*Had eye(nx) zeros(nx) ;
%         z*eye(nx)-Fad+K_k_ss*Had zeros(nx) zeros(nx) ;
%         N_k_ss*z -eye(nx) eye(nx)-N_k_ss*z
%         ];
%
%     H_yx(:,:,k)=Mat_L\Mat_R;
%
% end


% Hs=H_yx((nx*2+1):end,:,:);
