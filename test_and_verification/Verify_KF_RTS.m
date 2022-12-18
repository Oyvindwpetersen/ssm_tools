%% Verification of KF/RTS

clc
clear all
close all

%% Create system

mod=ImportSimplySupportQuick(2);

mod.dt=0.05;

%% Disc force

mod.a_cell={'10_U' '20_U' '30_U' '40_U' '50_U' '60_U' '70_U' '80_U' '90_U'}
mod.d_cell=mod.a_cell;
mod.p_cell={'30_U'}
[mod.Sd,mod.Sa,mod.Sp]=DofSelection(mod.d_cell,mod.a_cell,mod.p_cell,mod.doflabel);

[mod.A mod.B mod.G mod.J mod.Ac mod.Bc]=ssmod_modal(mod.phi,mod.Omega,mod.Gamma,mod.Sa,mod.Sd,mod.Sp,mod.dt,'force','disc');

mod.ny=size(mod.Sa,1); mod.nx=size(mod.A,1); mod.np=size(mod.Sp,2);

%%

mod.Q=eye(mod.nx)*[1e-2]^2;
mod.R=eye(mod.ny)*[1e-4]^2;
mod.S=zeros(mod.nx,mod.ny);

mod.S=[ones(mod.nx/2,mod.ny)*0.1 ; ones(mod.nx/2,mod.ny)*-0.1].*diag(mod.Q).^0.5.*(diag(mod.R).').^0.5;

Ctot=[mod.Q mod.S ; mod.S.' mod.R]
Ctot=plotcorr([mod.Q mod.S ; mod.S.' mod.R]);

sim.nt=1e6;
sim.t=[1:sim.nt]*mod.dt;

%% Generate noise

[sim.w,sim.v]=cov_noisegen(mod.Q,mod.R,mod.S,sim.t);

close all

sim.x0=zeros(mod.nx,1);
[sim.x,sim.y]=ssmod_forward_stoch(mod.A,[],mod.G,[],[],sim.x0,[],sim.w,sim.v);

plotTime(sim.t,sim.x);
plotTime(sim.t,sim.y);

tilefigs

%% Kalman filter and RTS smoother

P_0_0=[]

[x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KalmanFilter(mod.A,mod.G,mod.Q,mod.R,mod.S,sim.y,sim.x0,P_0_0,'forcescaling','yes');
[x_k_kb,x_k_kminb,P_k_kb,P_k_kminb,K_k_ssb]=KalmanFilter(mod.A,mod.G,mod.Q,mod.R,mod.S*0,sim.y,sim.x0,P_0_0,'forcescaling','yes');

% For RTS, S must be handled by decoupling and new system matrices

mod.A_star=mod.A-mod.S/mod.R*mod.G;
[x_k_N,P_k_N]=RTSSmoother(mod.A_star,x_k_k,x_k_kmin,P_k_k,P_k_kmin);
[x_k_Nb,P_k_Nb]=RTSSmoother(mod.A,x_k_kb,x_k_kminb,P_k_kb,P_k_kminb);
[x_k_N_false,P_k_N_false]=RTSSmoother(mod.A,x_k_k,x_k_kmin,P_k_k,P_k_kmin);

% close all
% plotTime(t,x,x_k_kmin,x_k_kminb);
% plotTime(t,x,x_k_k,x_k_kb);
% 
% plotTime(t,x,x_k_N,x_k_Nb,x_k_N_false);
% plotFreq(t,x-x_k_N,x-x_k_Nb,x-x_k_N_false,'xlim',[0 5]);

%%

clc
close all

var_filt=std(sim.x-x_k_k,0,2).^2
var_filt_b=std(sim.x-x_k_kb,0,2).^2

var_smooth=std(sim.x-x_k_N,0,2).^2
var_smooth_b=std(sim.x-x_k_Nb,0,2).^2
var_smooth_false=std(sim.x-x_k_N_false,0,2).^2

figure(); hold on; grid on; ylog;
title('Empirical state error variance'); legend show; 
plot(var_filt,'db','Markersize',8,'DisplayName','KF');
plot(var_filt_b,'*r','Markersize',8,'DisplayName','KF (S=0)');

plot(var_smooth,'om','Markersize',8,'DisplayName','RTS');
plot(var_smooth_b,'sk','Markersize',8,'DisplayName','RTS (S=0)');
plot(var_smooth_false,'xg','Markersize',8,'DisplayName','RTS (no correction)');

[P_k_kmin_ss,~,~,info]=idare(mod.A.',mod.G.',mod.Q,mod.R,mod.S); %,'noscaling'
P_k_k_ss=P_k_kmin_ss-P_k_kmin_ss*mod.G.'/(mod.G*P_k_kmin_ss*mod.G.'+mod.R)*mod.G*P_k_kmin_ss;

figure(); hold on; grid on; ylog;
title('P_{k|k-1} state error variance'); legend show; 
plot(diag(P_k_kmin),'db','Markersize',8,'DisplayName','KF');
plot(diag(P_k_kminb),'*r','Markersize',8,'DisplayName','KF (S=0)');
plot(diag(P_k_kmin_ss),'xg','Markersize',8,'DisplayName','DARE');

figure(); hold on; grid on; ylog;
title('P_{k|k} state error variance'); legend show; 
plot(diag(P_k_k),'db','Markersize',8,'DisplayName','KF');
plot(diag(P_k_kb),'*r','Markersize',8,'DisplayName','KF (S=0)');
plot(diag(P_k_k_ss),'xg','Markersize',8,'DisplayName','DARE');

figure(); hold on; grid on; ylog;
title('Ratio state error variance empirical/predicted'); legend show; 
plot(var_filt./diag(P_k_k),'db','Markersize',8,'DisplayName','KF');
plot(var_filt_b./diag(P_k_kb),'*r','Markersize',8,'DisplayName','KF (S=0)');

% figure(); hold on; grid on; ylog;
% title('Ratio state error variance empirical/predicted'); legend show; 
plot(var_smooth./diag(P_k_N),'om','Markersize',8,'DisplayName','RTS');
plot(var_smooth_b./diag(P_k_Nb),'sk','Markersize',8,'DisplayName','RTS (S=0)');
plot(var_smooth_false./diag(P_k_N_false),'xg','Markersize',8,'DisplayName','RTS (no correction)');

tilefigs([3 3],'l');
% return

return

%% Check JIS and smoother

close all

rng(1);

sim.p=BandLimitedWhiteNoise(0.1,3,sim.t,mod.np)*1;
sim.p=BellShapedNoise(0.1,1,sim.t,mod.np)*1;

[sim.w,sim.v]=cov_noisegen(mod.Q,mod.R,mod.S,sim.t);

[sim.x,sim.y]=ssmod_forward_stoch(mod.A,mod.B,mod.G,mod.J,[],sim.x0,sim.p,sim.w,sim.v);

plotTime(sim.t,sim.p); xlim([0 100]);
plotTime(sim.t,sim.x);
plotTime(sim.t,sim.y);

plotFreq(sim.t,sim.y);

clc
close all

P_0_1=eye(size(mod.Q));

[x_jis p_jis Px_jis_ss Pp_jis_ss ] = JIS_trunc_ss(mod.A,mod.B,mod.G,mod.J,sim.y,sim.x0,mod.R,mod.Q,mod.S,P_0_1,'trunc','no');

var_state_jis=std(sim.x-x_jis,0,2).^2
var_force_jis=std(sim.p-p_jis,0,2).^2

figure(); hold on; grid on; %ylog;
title('Ratio variance empirical/predicted'); legend show; 
plot(var_state_jis./diag(Px_jis_ss),'db','Markersize',8,'DisplayName','JIS (state)');
plot(var_force_jis./diag(Pp_jis_ss),'*r','Markersize',8,'DisplayName','JIS (force)');


L=10
[x_smooth p_smooth Px_smooth_ss Pp_smooth_ss ] = JIS_smooth(mod.A,mod.B,mod.G,mod.J,sim.y,mod.R,mod.Q,mod.S,sim.x0,P_0_1,L,'convtol',1e-8);

var_state_smooth=std(sim.x(:,1:end-L)-x_smooth(:,1:end-L),0,2).^2
var_force_smooth=std(sim.p(:,1:end-L)-p_smooth(:,1:end-L),0,2).^2


figure(); hold on; grid on; %ylog;
title('Ratio variance empirical/predicted'); legend show; 
plot(var_state_smooth./diag(Px_smooth_ss),'db','Markersize',8,'DisplayName','Smoother (state)');
plot(var_force_smooth./diag(Pp_smooth_ss),'*r','Markersize',8,'DisplayName','Smoother (force)');

tilefigs([3 3],'l');

plotTime(sim.t,sim.x,x_jis,x_smooth); xlimall([1000 2000]);
plotTime(sim.t,sim.p,p_jis,p_smooth); xlimall([1000 2000]);

tilefigs
% return

%% Check LFM

clc
close all

mod.sigma_p=std(sim.p(1,:),0,2);
mod.lambda=0.1;

mod.Gc=mod.G;
mod.Jc=mod.J;

[mod.Fc,mod.Lc,mod.Hc,mod.sigma_w]=ssmod_matern(mod.lambda,0,mod.sigma_p);

% mod.Fc=blockDiagonal(repcell(Fc,1,np));
% mod.Lc=blockDiagonal(repcell(Lc,1,np));
% mod.Hc=blockDiagonal(repcell(Hc,1,np));

mod.Qxd=blkdiag(mod.Q,zeros(size(mod.Fc,1)));

[mod.Fac,mod.Bac,mod.Hac,mod.Jac,mod.Fad,mod.Bad,mod.Had,mod.Jad,mod.Qad]=ssmod_lfm_aug(mod.Ac,mod.Bc,mod.Gc,mod.Jc,mod.Fc,mod.Hc,mod.Lc,mod.Qxd,mod.sigma_w,mod.dt);

% x0=zeros(nm*2,1);
% [x,y]=ssmod_forward_stoch(A,B,G,J,[],x0,p,w,v);
[sim.x,sim.y]=ssmod_forward_stoch(mod.A,mod.B,mod.G,mod.J,[],sim.x0,sim.p,sim.w,sim.v);

P_0_0=[]
sim.xa0=zeros(mod.nm*2+size(mod.Fc,1),1);
mod.S2=zeros(mod.nm*2+size(mod.Fc,1),mod.ny);
[xa_k_k,xa_k_kmin,Pa_k_k,Pa_k_kmin]=KalmanFilter(mod.Fad,mod.Had,mod.Qad,mod.R,mod.S2,sim.y,sim.xa0,P_0_0,'forcescaling','yes');

x_k_k=xa_k_k(1:mod.nm*2,:);
p_k_k=mod.Hc*xa_k_k((mod.nm*2+1):end,:);

[xa_k_N,Pa_k_N]=RTSSmoother(mod.Fad,xa_k_k,xa_k_kmin,Pa_k_k,Pa_k_kmin);

x_k_N=xa_k_N(1:mod.nm*2,:);
p_k_N=mod.Hc*xa_k_N((mod.nm*2+1):end,:);

plotTime(sim.t,sim.x,x_k_k,x_k_N);
plotTime(sim.t,sim.p,p_k_k,p_k_N);
plotFreq(sim.t,sim.p,p_k_k,p_k_N);


var_state_kf=std(sim.x-x_k_k,0,2).^2
var_force_kf=std(sim.p-p_k_k,0,2).^2

var_state_rts=std(sim.x-x_k_N,0,2).^2
var_force_rts=std(sim.p-p_k_N,0,2).^2

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

mod.sigma_p=std(sim.p(1,:),0,2);
mod.lambda=0.1;

lambda_vec=logspace(-3,2,100);
sigma_p_vec=logspace(-2,3,100);

for ind=1:100

% mod.lambda=lambda_vec(ind);
mod.sigma_p=sigma_p_vec(ind);

[mod.Fc,mod.Lc,mod.Hc,mod.sigma_w]=ssmod_matern(mod.lambda,0,mod.sigma_p);

[mod.Fac,mod.Bac,mod.Hac,mod.Jac,mod.Fad,mod.Bad,mod.Had,mod.Jad,mod.Qad]=ssmod_lfm_aug(mod.Ac,mod.Bc,mod.Gc,mod.Jc,mod.Fc,mod.Hc,mod.Lc,mod.Qxd,mod.sigma_w,mod.dt);

[xa_k_k,xa_k_kmin,Pa_k_k,Pa_k_kmin]=KalmanFilter(mod.Fad,mod.Had,mod.Qad,mod.R,mod.S2,sim.y,sim.xa0,P_0_0,'forcescaling','yes');
x_k_k=xa_k_k(1:mod.nm*2,:);
p_k_k=mod.Hc*xa_k_k((mod.nm*2+1):end,:);

[xa_k_N,Pa_k_N]=RTSSmoother(mod.Fad,xa_k_k,xa_k_kmin,Pa_k_k,Pa_k_kmin);
x_k_N=xa_k_N(1:mod.nm*2,:);
p_k_N=mod.Hc*xa_k_N((mod.nm*2+1):end,:);

var_state_kf=std(sim.x-x_k_k,0,2).^2;
var_force_kf=std(sim.p-p_k_k,0,2).^2;

var_state_rts=std(sim.x-x_k_N,0,2).^2;
var_force_rts=std(sim.p-p_k_N,0,2).^2;


Omega=mod.Had*Pa_k_kmin*mod.Had.'+mod.R;
[negLL(ind),negLL_1(ind),negLL_2(ind)]=loglik_evidence(mod.Had,xa_k_kmin,sim.y,Omega,'cut',[100 100]);

Vx_kf(ind,:)=var_state_kf;
Vp_kf(ind,:)=var_force_kf;

Vx_rts(ind,:)=var_state_rts;
Vp_rts(ind,:)=var_force_rts;

end


close all

figure(); hold on; grid on; ylog; xlog;
title('Variance KF state');
plot(lambda_vec,Vx_kf,'-','DisplayName','KF');
plot(lambda_vec,Vx_rts,'--','DisplayName','RTS');
legend show

figure(); hold on; grid on; ylog; xlog;
title('Variance KF force');
plot(lambda_vec,Vp_kf,'-','DisplayName','KF');
plot(lambda_vec,Vp_rts,'--','DisplayName','RTS');
legend show

figure(); hold on; grid on; xlog;
title('LogL');
plot(lambda_vec,negLL,'-','DisplayName','Sum');
plot(lambda_vec,negLL_1,'--','DisplayName','Model complexity');
plot(lambda_vec,negLL_2,'*-','DisplayName','Data fit');
legend show

tilefigs([3 3],'l');

