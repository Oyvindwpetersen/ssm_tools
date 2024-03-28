%%

clc
clear all
close all

%% Create system

mod=importbeamquick(6);

mod.dt=0.01;
% mod.dt=0.025/10;

%% Disc force

% mod.a_cell={'10_U' '40_U'  '80_U'  }
% mod.d_cell={'15_U'}

% mod.a_cell=genlabels('U',[10:10:90]); % For test modal expansion

mod.a_cell={'10_U'  '50_U'  '80_U' }
mod.d_cell={'15_U'  '55_U'}
mod.p_cell={'30_U'}
mod.e_cell={'25_U' '45_U' '80_U'}

[mod.Sd,mod.Sa,mod.Sp]=DofSelection(mod.d_cell,mod.a_cell,mod.p_cell,mod.doflabel);
[~,~,mod.Se]=DofSelection({},{},mod.e_cell,mod.doflabel);

[mod.A mod.B mod.G mod.J mod.Ac mod.Bc mod.Gc mod.Jc]=ssmod_modal(mod.phi,mod.Omega,mod.Gamma,mod.Sa,mod.Sd,mod.Sp,mod.dt,'force','disc');
[~, mod.Be ,~,mod.Je]=ssmod_modal(mod.phi,mod.Omega,mod.Gamma,mod.Sa,mod.Sd,mod.Se,mod.dt,'force','disc');

mod.ny=size(mod.Sa,1); mod.nx=size(mod.A,1); mod.np=size(mod.Sp,2);

%%

Ce=eye(3);
mod.R0=eye(mod.ny)*[1e-4]^2;
mod.Q0=eye(mod.nx)*[1e-8]^2;

mod.Q=mod.Be*Ce*mod.Be.'+mod.Q0;
mod.R=mod.Je*Ce*mod.Je.'+mod.R0;
mod.S=mod.Be*Ce*mod.Je.';

mod.Q=mod.Q*100;
mod.R=mod.R*100;
mod.S=mod.S*100;

plotcorr([mod.Q mod.S ; mod.S.' mod.R]);

sim.nt=1e5+1;
sim.t=[1:sim.nt]*mod.dt;

%% Generate noise

close all
% sim.p=BandLimitedWhiteNoise(0.1,3,sim.t,mod.np)*1;
sim.p=bellshapednoise(0.05,1,sim.t,mod.np)*10;

[sim.w,sim.v]=cov_noisegen(mod.Q,mod.R,mod.S,sim.t);

sim.x0=zeros(mod.nx,1);
[sim.x,sim.y]=ssmod_forward_stoch(mod.A,mod.B,mod.G,mod.J,[],sim.x0,sim.p,sim.w,sim.v);

% plottime(sim.t,sim.x);
plottime(sim.t,sim.y);

tilefigs

return

%% KF

[x_k_k,x_k_kmin,P_k_k,P_k_kmin]=KF(mod.A,mod.B,mod.G,mod.J,mod.Q,mod.R,mod.S,sim.y,sim.p,[],[],'steadystate',true);

close all

plottime(sim.t,sim.x,x_k_kmin,x_k_k);
plotfreq(sim.t,sim.x,x_k_kmin,x_k_k,'xlim',[0 10]);

% plotstat(diag(P_k_kmin),diag(P_k_k),std(sim.x-x_k_kmin,0,2).^2,std(sim.x-x_k_k,0,2).^2);
plotstat(std(sim.x-x_k_kmin,0,2).^2./diag(P_k_kmin),std(sim.x-x_k_k,0,2).^2./diag(P_k_k));


%% JIS

close all

mod.P01=[];

[x_jis_old p_jis_old Px_jis_ss_old Pp_jis_ss_old ] = JIS_ss_old(mod.A,mod.B,mod.G,mod.J,sim.y,sim.x0,mod.Q,mod.R,mod.S,mod.P01,'trunc',false);

[x_jis p_jis Px_jis_ss Pp_jis_ss] = JIS_ss(mod.A,mod.B,mod.G,mod.J,sim.y,sim.x0,mod.Q,mod.R,mod.S,mod.P01,'trunc',false);

fid.L=10;
[x_smooth p_smooth Px_smooth_ss Pp_smooth_ss ]=JIS_smooth(mod.A,mod.B,mod.G,mod.J,sim.y,mod.Q,mod.R,mod.S,sim.x0,mod.P01,fid.L,'convtol',1e-10);

plottime(sim.t,sim.p,p_jis_old,p_jis,p_smooth,'displayname',{'True' 'Filt old' 'Filt' 'Smooth'});

% close all

plotopt=struct();
plotopt.gap=[0.1 0.1];
plotopt.marg_h=[0.1 0.25];
plotopt.ylabel='Error SD';
plotopt.markersize=10;
plotopt.color={[0 0 1] [1 0 0]   [0 0 1] [1 0 0]  };
% plotopt.split={};

% States

plot_std_emp_x=sdstatemp(sim.x,{[x_jis] , [x_smooth]},[50 50])

plot_std_pred_x={};
plot_std_pred_x{1}=[diag(Px_jis_ss).^0.5 ];
plot_std_pred_x{2}=[diag(Px_smooth_ss).^0.5 ];

plotstat(plot_std_emp_x{1}./plot_std_pred_x{1},plot_std_emp_x{2}./plot_std_pred_x{2},plotopt);

% legendadjust('east',[-0.1 0.05],4,12);
% axistight(gcf,[0.3 1],'ylog','keepx');

% Force
plot_std_emp_p=sdstatemp(sim.p,{[p_jis] , [p_smooth]},[50 50])

plot_std_pred_p{1}=[diag(Pp_jis_ss).^0.5];
plot_std_pred_p{2}=[diag(Pp_smooth_ss).^0.5];

plotstat(plot_std_emp_p{1},plot_std_pred_p{1},plotopt);

plotstat(plot_std_emp_p{2},plot_std_pred_p{2},plotopt);

tilefigs

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


%% KF

[x_k_k,x_k_kmin,x_k_N,P_k_k,P_k_kmin,P_k_N]=KF_RTS(mod.A,mod.B,mod.G,mod.J,mod.Q,mod.R,mod.S,sim.y,sim.p);

% Check statistics (best for long time series, band limited white noise input)

close all

plotopt=struct();
plotopt.gap=[0.1 0.1];
plotopt.marg_h=[0.1 0.25];
plotopt.ylabel='Error SD';
plotopt.displayname={'Pred' 'Filt' 'Smooth' };
plotopt.markersize=10;
plotopt.color={[0 0 1] [1 0 0] [0 0 0]  };
plotopt.marker={ 'x' 'o' 'd'};

% States

plot_std_emp_x=sdstatemp(sim.x,{x_k_kmin , x_k_k , x_k_N},[50 50])

plot_std_pred_x={};
plot_std_pred_x{1}=[diag(P_k_kmin).^0.5];
plot_std_pred_x{2}=[diag(P_k_k).^0.5];
plot_std_pred_x{3}=[diag(P_k_N).^0.5];

plotstat(plot_std_emp_x{1},plot_std_pred_x{1},plotopt);
plotstat(plot_std_emp_x{2},plot_std_pred_x{2},plotopt);
plotstat(plot_std_emp_x{3},plot_std_pred_x{3},plotopt);

plotopt.normalize=true;
plotstat(plot_std_emp_x{:},plotopt);

tilefigs
