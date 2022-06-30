function [Hp Hf Hs]=KF_RTS_SS_TF(Fad,Had,Q,R,S,dt,w_axis)

%%

%% Run KF and RTS

x0=zeros(size(Fad,1),1);
P_0_0=[];

y=zeros(size(Had,1),1e3);

% Obtain matrices
[x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KalmanFilter(Fad,Had,Q,R,S,y,x0,P_0_0,'steadystate','yes');
[x_k_N,P_k_N]=RTSSmoother(Fad,x_k_k,x_k_kmin,P_k_k,P_k_kmin,'steadystate','yes');

if any(any(S))
    error('Implementation not supporting S~=0, due to modification necessary in RTS smoother');
end

%% TF of filter and smoother

P_k_kmin_ss=P_k_kmin;
P_k_k_ss=P_k_k;

Omega_k_ss=(Had*P_k_kmin_ss*Had.'+R); Omega_k_ss=forcesym(Omega_k_ss);
Omega_k_ss_inv=eye(size(Omega_k_ss))/Omega_k_ss; Omega_k_ss_inv=forcesym(Omega_k_ss_inv);
K_k_ss=(Fad*P_k_kmin_ss*Had.'+S)*Omega_k_ss_inv;

N_k_ss=P_k_k_ss*Fad.'/P_k_kmin_ss;

nx=size(Fad,1);
ny=size(Had,1);

H_yx=zeros(nx*3,ny,length(w_axis));

Mat_R=[P_k_kmin_ss*Had.'*Omega_k_ss_inv ; K_k_ss ; zeros(nx,ny) ];

for k=1:length(w_axis)
    z=exp(1i*w_axis(k)*dt);
    
%     Mat=[
%         -eye(nx)+K_k_ss*Had eye(nx) zeros(nx) ;
%         z*eye(nx) -Fad zeros(nx) ;
%         N_k_ss*z -eye(nx) eye(nx)-N_k_ss*z 
%         ];
    
    Mat_L=[
        -eye(nx)+P_k_kmin_ss*Had.'*Omega_k_ss_inv*Had eye(nx) zeros(nx) ;
        z*eye(nx)-Fad+K_k_ss*Had zeros(nx) zeros(nx) ;
        N_k_ss*z -eye(nx) eye(nx)-N_k_ss*z
        ];
    
    H_yx(:,:,k)=Mat_L\Mat_R;

end

Hp=H_yx(1:nx,:,:);
Hf=H_yx((nx+1):(nx*2),:,:);
Hs=H_yx((nx*2+1):end,:,:);
