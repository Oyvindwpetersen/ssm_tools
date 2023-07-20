%%

clc
clear all
close all

%% Create system

mod=importbeamquick(6);

mod.dt=0.05/10;

%% Disc force

mod.a_cell={'10_U'  '50_U'   '80_U'}
mod.d_cell={'15_U'  '55_U'}
mod.p_cell={'30_U'}
mod.e_cell={'25_U' '45_U' '65_U'}

[mod.Sd,mod.Sa,mod.Sp]=DofSelection(mod.d_cell,mod.a_cell,mod.p_cell,mod.doflabel);
[~,~,mod.Se]=DofSelection({},{},mod.e_cell,mod.doflabel);

[mod.A mod.B mod.G mod.J mod.Ac mod.Bc]=ssmod_modal(mod.phi,mod.Omega,mod.Gamma,mod.Sa,mod.Sd,mod.Sp,mod.dt,'force','disc');
[~, mod.Be ,~,mod.Je]=ssmod_modal(mod.phi,mod.Omega,mod.Gamma,mod.Sa,mod.Sd,mod.Se,mod.dt,'force','disc');

mod.ny=size(mod.Sa,1); mod.nx=size(mod.A,1); mod.np=size(mod.Sp,2);

%%

Ce=eye(3);
mod.R0=eye(mod.ny)*[1e-4]^2;
mod.Q0=eye(mod.nx)*[1e-6]^2;

mod.Q=mod.Be*Ce*mod.Be.'+mod.Q0;
mod.R=mod.Je*Ce*mod.Je.'+mod.R0;
mod.S=mod.Be*Ce*mod.Je.';
% Ctot=[mod.Q mod.S ; mod.S.' mod.R]
Ccorr=plotcorr([mod.Q mod.S ; mod.S.' mod.R]);

sim.nt=1e5;
sim.t=[1:sim.nt]*mod.dt;

%% Generate noise

close all
% sim.p=BandLimitedWhiteNoise(0.1,3,sim.t,mod.np)*1;
sim.p=BellShapedNoise(0.1,1,sim.t,mod.np)*10;

[sim.w,sim.v]=cov_noisegen(mod.Q,mod.R,mod.S,sim.t);

sim.x0=zeros(mod.nx,1);
[sim.x,sim.y]=ssmod_forward_stoch(mod.A,mod.B,mod.G,mod.J,[],sim.x0,sim.p,sim.w,sim.v);

% plottime(sim.t,sim.x);
plottime(sim.t,sim.y);

tilefigs


%% Estimate

close all

% mod.P01=eye(size(mod.Q));
mod.P01=[];

[x_jis p_jis Px_jis_ss Pp_jis_ss ] = JIS_ss(mod.A,mod.B,mod.G,mod.J,sim.y,sim.x0,mod.Q,mod.R,mod.S,mod.P01,'trunc',false);

fid.L=10;
[x_smooth p_smooth Px_smooth_ss Pp_smooth_ss ]=JIS_smooth(mod.A,mod.B,mod.G,mod.J,sim.y,mod.Q,mod.R,mod.S,sim.x0,mod.P01,fid.L,'convtol',1e-10);

plottime(sim.t,sim.p,p_jis,p_smooth,'displayname',{'True' 'Filt' 'Smooth'});

plotfreq(sim.t,sim.p,p_jis,p_smooth,'xlim',[0 50],'displayname',{'True' 'Filt' 'Smooth'});

% plottime(sim.t,sim.x,x_jis,x_smooth);

% plotfreq(sim.t,sim.x(1:6,:),x_jis(1:6,:),x_smooth(1:6,:),'xlim',[0 50]);

tilefigs


% return
%% Check statistics (best for long time series, band limited white noise input)

close all

plotopt=struct();
plotopt.gap=[0.1 0.1];
plotopt.marg_h=[0.1 0.25];
plotopt.ylabel='Error SD';
plotopt.displayname={'JIS' 'Smooth' };
plotopt.displayname2={'JIS' 'Smooth' };
plotopt.markersize=10;
plotopt.color={[0 0 1] [1 0 0]   [0 0 1] [1 0 0]  };
% plotopt.marker={ 'x' 'x' 'x' 'x'  };
% plotopt.xticklabel={'x_1' 'x_2' 'x_3' 'x_4' 'x_5' 'x_6'  'x_7' 'x_8' 'x_9' 'x_{10}' 'p'};
% plotopt.split={};

% States

plot_std_emp_x=sdstatemp(sim.x,{[x_jis] , [x_smooth]},[50 50])

plot_std_pred_x={};
plot_std_pred_x{1}=[diag(Px_jis_ss).^0.5 ];
plot_std_pred_x{2}=[diag(Px_smooth_ss).^0.5 ];

plotstatstem(plot_std_emp_x,plot_std_pred_x,plotopt);

legendadjust('east',[-0.1 0.05],4,12);
axistight(gcf,[0.3 1],'ylog','keepx');

% Force
plot_std_emp_p=sdstatemp(sim.p,{[p_jis] , [p_smooth]},[50 50])

plot_std_pred_p={};
plot_std_pred_p{1}=[diag(Pp_jis_ss).^0.5];
plot_std_pred_p{2}=[diag(Pp_smooth_ss).^0.5];

plotstatstem(plot_std_emp_p,plot_std_pred_p,plotopt);

return

%% Test lag

% L_mat=[1 3 5 10 20 30 40 50 60];
L_mat=[10:10:30];

for k=1:length(L_mat)
    [~,~,P_x_ss,P_p_ss ]=JIS_smooth(mod.A,mod.B,mod.G,mod.J,sim.y,mod.Q,mod.R,mod.S,sim.x0,mod.P01,L_mat(k));
    Px_all(k,:)=diag(P_x_ss);
    Pp_all(k,:)=diag(P_p_ss);
end

close all
figure(); hold on; grid on;
plot(L_mat,Px_all); ylog;

figure(); hold on; grid on;
plot(L_mat,Pp_all); ylog;


%%

P_0_0=[]

[x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KalmanFilterWithInput(mod.A,mod.B,mod.G,mod.J,mod.Q,mod.R,mod.S,sim.y,sim.p,sim.x0,P_0_0);


mod.A_star=mod.A-mod.S/mod.R*mod.G;
mod.B_star=[mod.B-mod.S/mod.R*mod.J mod.S/mod.R];
mod.G_star=[mod.G ];
mod.J_star=[mod.J zeros(mod.ny,mod.ny)];

sim.p_star=[sim.p ; sim.y];
mod.Q_star=mod.Q-mod.S/mod.R*mod.S.';
mod.S_star=zeros(size(mod.A,1),size(mod.G,1));
mod.R_star=mod.R;

[x_k_k_tr,x_k_kmin_tr,P_k_k_tr,P_k_kmin_tr]=KalmanFilterWithInput(mod.A_star,mod.B_star,mod.G_star,mod.J_star,mod.Q_star,mod.R_star,mod.S_star,sim.y,sim.p_star,sim.x0,P_0_0);

close all
plottime(sim.t,x_k_k,x_k_k_tr)
plottime(sim.t,x_k_kmin,x_k_kmin_tr)

[x_k_N_tr,P_k_N_tr]=RTSSmoother(mod.A_star,x_k_k_tr,x_k_kmin_tr,P_k_k_tr,P_k_kmin_tr);

% For RTS, S must be handled by decoupling and new system matrices
% [x_k_kb,x_k_kminb,P_k_kb,P_k_kminb,K_k_ssb]=KalmanFilter(mod.A,mod.G,mod.Q,mod.R,mod.S*0,sim.y,sim.x0,P_0_0,'forcescaling','yes');

% mod.A_star=mod.A-mod.S/mod.R*mod.G;
% [x_k_N,P_k_N]=RTSSmoother(mod.A_star,x_k_k,x_k_kmin,P_k_k,P_k_kmin);
% [x_k_Nb,P_k_Nb]=RTSSmoother(mod.A,x_k_kb,x_k_kminb,P_k_kb,P_k_kminb);
% [x_k_N_false,P_k_N_false]=RTSSmoother(mod.A,x_k_k,x_k_kmin,P_k_k,P_k_kmin);

%% Check statistics (best for long time series, band limited white noise input)

close all

plotopt=struct();
plotopt.gap=[0.1 0.1];
plotopt.marg_h=[0.1 0.25];
plotopt.ylabel='Error SD';
plotopt.displayname={'Pred' 'Filt' 'Smooth' };
plotopt.displayname2={'Pred' 'Filt' 'Smooth' };
plotopt.markersize=10;
plotopt.color={[0 0 1] [1 0 0] [0 0 0]  [0 0 1] [1 0 0]  [0 0 0] };
plotopt.marker={ 'x' 'x' 'x' 'x'  };
% plotopt.split={};

% States

plot_std_emp_x=sdstatemp(sim.x,{x_k_kmin_tr , x_k_k_tr , x_k_N_tr},[50 50])

plot_std_pred_x={};
plot_std_pred_x{1}=[diag(P_k_kmin_tr).^0.5 ];
plot_std_pred_x{2}=[diag(P_k_k_tr).^0.5 ];
plot_std_pred_x{3}=[diag(P_k_N_tr).^0.5 ];

plotstatstem(plot_std_emp_x,plot_std_pred_x,plotopt);

legendadjust('east',[-0.1 0.05],4,12);
axistight(gcf,[0.3 1],'ylog','keepx');

return
