%% Verification of KF/RTS


clc
clear all
close all

% rng('Default');

%% Create system

nm=4;
dt=0.05;

Omega=diag([0.374 1.44 2.5 3.636]*2*pi);
Xi=diag([1.86 1.38 0.92 1.2]/100);
Gamma=2*Omega*Xi;

Omega=Omega(1:nm,1:nm);
Xi=Xi(1:nm,1:nm);
Gamma=Gamma(1:nm,1:nm);

% Mode shapes
L=100
x_axis=linspace(0,L,100);
for k=1:nm
    phi(:,k)=sin(k*pi*x_axis/L)/1e4;
end

doflabel=getLabel('U1',[1:length(phi)]);

%%% Modal
% a_cell={'20_U1' '30_U1' '40_U1' '60_U1'}
% d_cell=a_cell;
% p_cell={'20_U1' '60_U1' '80_U1'     '10_U1' '30_U1' '40_U1' '50_U1' '90_U1'}
% [Sd,Sa,Sp]=DofSelection(d_cell,a_cell,p_cell,doflabel);
% 
% [A B G J Ac Bc]=ssmod_modal(phi,Omega,Gamma,Sa,Sd,Sp,dt,'force','disc');
% 
% ny=size(Sa,1); nx=size(A,1); np=nm;
% 
% phi_p=Sp.'*phi;
% Cf=1e8*eye(size(phi_p,1));
% Cp=phi_p.'*Cf*phi_p
% plotcorr(Cp)
% 
% R0=eye(ny)*[1e-4]^2;
% Q=B*Cp*B.'; Q=Q+eye(size(Q))*max(max(Q))*0.01;
% R=J*Cp*J.'+R0
% S=B*Cp*J.'*0

%%% Disc
a_cell={'10_U1' '20_U1' '30_U1' '40_U1' '50_U1' '60_U1' '70_U1' '80_U1' '90_U1'}
d_cell=a_cell;
p_cell={'20_U1'}% '60_U1' '80_U1'  '90_U1'
[Sd,Sa,Sp]=DofSelection(d_cell,a_cell,p_cell,doflabel);

[A B G J Ac Bc]=ssmod_modal(phi,Omega,Gamma,Sa,Sd,Sp,dt,'force','disc');

ny=size(Sa,1); nx=size(A,1); np=size(Sp,2);

R0=eye(ny)*[1e-8]^2;
Q=eye(size(A))*1e-12;
R=R0;
S=zeros(nx,ny);

Ctot=[Q S ; S.' R]
Ccorr=plotcorr(Ctot);

nt=1e4;
t=[1:nt]*dt;

%% Generate noise

[w,v]=cov_noisegen(Q,R,S,t);

%%

close all

x0=zeros(nx,1);
[x,y]=ssmod_forward_stoch(A,[],G,[],[],x0,[],w,v);

% plotTime(t,x,x_true);
% plotTime(t,y,y_true);

tilefigs

%% Kalman filter and RTS smoother

P_0_0=[]

[x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KalmanFilter(A,G,Q,R,S,y,x0,P_0_0,'forcescaling','yes');
[x_k_kb,x_k_kminb,P_k_kb,P_k_kminb,K_k_ssb]=KalmanFilter(A,G,Q,R,S*0,y,x0,P_0_0,'forcescaling','yes');

% For RTS, S must be handled by decoupling and new system matrices

A_star=A-S/R*G;
[x_k_N,P_k_N]=RTSSmoother(A_star,x_k_k,x_k_kmin,P_k_k,P_k_kmin);
[x_k_Nb,P_k_Nb]=RTSSmoother(A_star,x_k_kb,x_k_kminb,P_k_kb,P_k_kminb);
[x_k_N_false,P_k_N_false]=RTSSmoother(A,x_k_k,x_k_kmin,P_k_k,P_k_kmin);

% close all
% plotTime(t,x,x_k_kmin,x_k_kminb);
% plotTime(t,x,x_k_k,x_k_kb);
% 
% plotTime(t,x,x_k_N,x_k_Nb,x_k_N_false);
% plotFreq(t,x-x_k_N,x-x_k_Nb,x-x_k_N_false,'xlim',[0 5]);

%%

clc
close all

var_filt=std(x-x_k_k,0,2).^2
var_filt_b=std(x-x_k_kb,0,2).^2

var_smooth=std(x-x_k_N,0,2).^2
var_smooth_b=std(x-x_k_Nb,0,2).^2
var_smooth_false=std(x-x_k_N_false,0,2).^2

figure(); hold on; grid on; ylog;
title('Empirical variance'); legend show; 
plot(var_filt,'db','Markersize',8,'DisplayName','KF');
plot(var_filt_b,'*r','Markersize',8,'DisplayName','KF (S=0)');

plot(var_smooth,'om','Markersize',8,'DisplayName','RTS');
plot(var_smooth_b,'sk','Markersize',8,'DisplayName','RTS (S=0)');
plot(var_smooth_false,'xg','Markersize',8,'DisplayName','RTS (no corr)');

[P_k_kmin_ss,~,~,info]=idare(A.',G.',Q,R,S); %,'noscaling'
P_k_k_ss=P_k_kmin_ss-P_k_kmin_ss*G.'/(G*P_k_kmin_ss*G.'+R)*G*P_k_kmin_ss;

figure(); hold on; grid on; ylog;
title('P_{k|k-1} variance'); legend show; 
plot(diag(P_k_kmin),'db','Markersize',8,'DisplayName','KF');
plot(diag(P_k_kminb),'*r','Markersize',8,'DisplayName','KF (S=0)');
plot(diag(P_k_kmin_ss),'xg','Markersize',8,'DisplayName','DARE');

figure(); hold on; grid on; ylog;
title('P_{k|k} variance'); legend show; 
plot(diag(P_k_k),'db','Markersize',8,'DisplayName','KF');
plot(diag(P_k_kb),'*r','Markersize',8,'DisplayName','KF (S=0)');
plot(diag(P_k_k_ss),'xg','Markersize',8,'DisplayName','DARE');

figure(); hold on; grid on; ylog;
title('Ratio variance empirical/predicted'); legend show; 
plot(var_filt./diag(P_k_k),'db','Markersize',8,'DisplayName','KF');
plot(var_filt_b./diag(P_k_kb),'*r','Markersize',8,'DisplayName','KF (S=0)');


figure(); hold on; grid on; ylog;
title('Ratio variance empirical/predicted'); legend show; 
plot(var_smooth./diag(P_k_N),'db','Markersize',8,'DisplayName','RTS');
plot(var_smooth_b./diag(P_k_Nb),'*r','Markersize',8,'DisplayName','RTS (S=0)');
plot(var_smooth_false./diag(P_k_N_false),'xg','Markersize',8,'DisplayName','RTS (no corr)');

tilefigs([3 3],'l');
% return

%% Check JIS and smoother

close all

rng(1);

p=BandLimitedWhiteNoise(0.3,4,t,np)*100;
% p=BandLimitedWhiteNoise(3,4,t,np)*50;

% p=trianglepulse(3,3.25,3.3,t)*1000;

[w,v]=cov_noisegen(Q,R,S,t);

w=w*0;

x0=zeros(nx,1);
[x,y]=ssmod_forward_stoch(A,B,G,J,[],x0,p,w,v);
% [x2,y2]=ssmod_forward_stoch(A,B,G,J,[],x0,p*0,w,v);

clc
close all

P_0_1=eye(size(Q));

[x_jis p_jis Px_jis_ss Pp_jis_ss ] = JIS_trunc_ss(A,B,G,J,y,x0,R,Q,S,P_0_1,'trunc','no');

var_state_jis=std(x-x_jis,0,2).^2
var_force_jis=std(p-p_jis,0,2).^2

figure(); hold on; grid on; %ylog;
title('Ratio variance empirical/predicted'); legend show; 
plot(var_state_jis./diag(Px_jis_ss),'db','Markersize',8,'DisplayName','JIS (state)');
plot(var_force_jis./diag(Pp_jis_ss),'*r','Markersize',8,'DisplayName','JIS (force)');

[x_smooth p_smooth Px_smooth_ss Pp_smooth_ss ] = JIS_smooth(A,B,G,J,y,R,Q,S,x0,P_0_1,10,'convtol',1e-8);

% var_state_smooth=std(x-x_smooth,0,2).^2
% var_force_smooth=std(p-p_smooth,0,2).^2
var_state_smooth=std(x(:,1:end-L)-x_smooth(:,1:end-L),0,2).^2
var_force_smooth=std(p(:,1:end-L)-p_smooth(:,1:end-L),0,2).^2


figure(); hold on; grid on; %ylog;
title('Ratio variance empirical/predicted'); legend show; 
plot(var_state_smooth./diag(Px_smooth_ss),'db','Markersize',8,'DisplayName','Smoother (state)');
plot(var_force_smooth./diag(Pp_smooth_ss),'*r','Markersize',8,'DisplayName','Smoother (force)');

tilefigs([3 3],'l');

plotTime(t,x,x_jis,x_smooth);
plotTime(t,p,p_jis,p_smooth);%xlimall([1000 2000]);
% 
% plotTime(t,y);


% return

%%

clc
close all

nt=1e4;
t=[1:nt]*dt;

p=BandLimitedWhiteNoise(0.5,0.7,t,np)*50;
% p=sin(2*pi*t)*50;

sigma_p=std(p(1,:),0,2);
lambda=0.1;

Gc=G;
Jc=J;

[Fc,Lc,Hc,sigma_w]=ssmod_matern(lambda,0,sigma_p);

Fc=blockDiagonal(repcell(Fc,1,np));
Lc=blockDiagonal(repcell(Lc,1,np));
Hc=blockDiagonal(repcell(Hc,1,np));

Qxd=eye(nm*2+size(Fc,1))*1e-6;

[Fac,Bac,Hac,Jac,Fad,Bad,Had,Jad,Qad]=ssmod_lfm_aug(Ac,Bc,Gc,Jc,Fc,Hc,Lc,Qxd,sigma_w,dt);

w=zeros(nm*2,nt);
[v]=mvnrnd(zeros(size(R0,1),1),R0,nt).';

x0=zeros(nm*2,1);
[x,y]=ssmod_forward_stoch(A,B,G,J,[],x0,p,w,v);

P_0_0=[]
xa0=zeros(nm*2+size(Fc,1),1);
S=zeros(nm*2+size(Fc,1),ny);
[xa_k_k,xa_k_kmin,Pa_k_k,Pa_k_kmin]=KalmanFilter(Fad,Had,Qad,R0,S,y,xa0,P_0_0,'forcescaling','yes');

x_k_k=xa_k_k(1:nm*2,:);
p_k_k=Hc*xa_k_k((2*nm+1):end,:);

[xa_k_N,Pa_k_N]=RTSSmoother(Fad,xa_k_k,xa_k_kmin,Pa_k_k,Pa_k_kmin);

x_k_N=xa_k_N(1:nm*2,:);
p_k_N=Hc*xa_k_N((2*nm+1):end,:);

plotTime(t,x,x_k_k,x_k_N);
plotTime(t,p,p_k_k,p_k_N);
plotFreq(t,p,p_k_k,p_k_N);


var_state_kf=std(x-x_k_k,0,2).^2
var_force_kf=std(p-p_k_k,0,2).^2

var_state_rts=std(x-x_k_N,0,2).^2
var_force_rts=std(p-p_k_N,0,2).^2

figure(); hold on; grid on; ylog;
title('Empirical variance'); legend show; 
plot(var_state_kf,'db','Markersize',8,'DisplayName','KF (state)');
plot(var_state_rts,'*r','Markersize',8,'DisplayName','RTS (state)');

figure(); hold on; grid on; ylog;
title('Empirical variance'); legend show; 
plot(var_force_kf,'db','Markersize',8,'DisplayName','KF (force)');
plot(var_force_rts,'*r','Markersize',8,'DisplayName','RTS (force)');


tilefigs([3 3],'l');

%%

lambda_vec=logspace(-3,2,100);

sigma_p_scale=logspace(-2,3,100);


for ind=1:length(lambda_vec)

% lambda=lambda_vec(ind);
lambda=0.1;

[Fc,Lc,Hc,sigma_w]=ssmod_matern(lambda,0,sigma_p*sigma_p_scale(ind));

Fc=blockDiagonal(repcell(Fc,1,np));
Lc=blockDiagonal(repcell(Lc,1,np));
Hc=blockDiagonal(repcell(Hc,1,np));

[Fac,Bac,Hac,Jac,Fad,Bad,Had,Jad,Qad]=ssmod_lfm_aug(Ac,Bc,Gc,Jc,Fc,Hc,Lc,Qxd,sigma_w,dt);

[xa_k_k,xa_k_kmin,Pa_k_k,Pa_k_kmin]=KalmanFilter(Fad,Had,Qad,R0,S,y,xa0,P_0_0,'forcescaling','yes');
x_k_k=xa_k_k(1:nm*2,:);
p_k_k=Hc*xa_k_k((2*nm+1):end,:);

[xa_k_N,Pa_k_N]=RTSSmoother(Fad,xa_k_k,xa_k_kmin,Pa_k_k,Pa_k_kmin);
x_k_N=xa_k_N(1:nm*2,:);
p_k_N=Hc*xa_k_N((2*nm+1):end,:);

var_state_kf=std(x-x_k_k,0,2).^2;
var_force_kf=std(p-p_k_k,0,2).^2;

var_state_rts=std(x-x_k_N,0,2).^2;
var_force_rts=std(p-p_k_N,0,2).^2;


Omega=Had*Pa_k_kmin*Had.'+R;
[negLL(ind),negLL_1(ind),negLL_2(ind)]=loglik_evidence(Had,xa_k_kmin,y,Omega,'cut',[100 100]);

Vx_kf(ind,:)=var_state_kf;
Vp_kf(ind,:)=var_force_kf;

Vx_rts(ind,:)=var_state_rts;
Vp_rts(ind,:)=var_force_rts;

end


close all

figure(); hold on; grid on; ylog; xlog;
title('Variance KF state');
plot(lambda_vec,Vx_kf,'-');
plot(lambda_vec,Vx_rts,'--');

figure(); hold on; grid on; ylog; xlog;
title('Variance KF force');
plot(lambda_vec,Vp_kf,'-');
plot(lambda_vec,Vp_rts,'--');


figure(); hold on; grid on; xlog;
title('LogL');
plot(lambda_vec,negLL,'-');
plot(lambda_vec,negLL_1,'--');
plot(lambda_vec,negLL_2,'*-');

tilefigs([3 3],'l');

