%%

clc
clear all
close all

%% Create system

mod=importbeamquick(8);

mod.dt=0.01;

%% Disc force

% mod.a_cell=genlabels('U',[10:10:90]); % For test modal expansion

mod.a_cell={'10_U'  '50_U'  } %'80_U' 
mod.d_cell={'15_U'  '55_U'}
mod.p_cell={'30_U'}
mod.e_cell={'25_U' '45_U' '80_U'}

% mod.a_cell={'50_U' }
% mod.d_cell={'15_U'  '55_U'}
% mod.p_cell={'30_U' '40_U'}
% mod.e_cell={'25_U' '45_U' '80_U'}

% mod.a_cell={'50_U' '15_U'  '55_U'}
% mod.d_cell={}
% mod.p_cell={'30_U' '40_U'}
% mod.e_cell={'25_U' '45_U' '80_U'}

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

% mod.Q=mod.Q/10000;
% mod.R=mod.R/10000;
% mod.S=mod.S/10000;

plotcorr([mod.Q mod.S ; mod.S.' mod.R]);

sim.nt=1e6+1;
sim.t=[1:sim.nt]*mod.dt;

%% Generate noise

close all
% sim.p=bellshapednoise(0.05,1,sim.t,mod.np)*10;

sim.p(1,:)=bellshapednoise(0.05,1,sim.t,1)*10;
% sim.p(2,:)=bellshapednoise(0.05,1,sim.t,1)*10;

[sim.w,sim.v]=cov_noisegen(mod.Q,mod.R,mod.S,sim.t);

sim.x0=zeros(mod.nx,1);
[sim.x,sim.y]=ssmod_forward_stoch(mod.A,mod.B,mod.G,mod.J,[],sim.x0,sim.p,sim.w,sim.v);

% plottime(sim.t,sim.x);
plottime(sim.t,sim.y);

tilefigs

%% JIS

close all

mod.P01=[];

[x_jis_old p_jis_old Px_jis_ss_old Pp_jis_ss_old] = JIS_ss_old(mod.A,mod.B,mod.G,mod.J,sim.y,sim.x0,mod.Q,mod.R,mod.S,mod.P01,'trunc',false);

[x_jis p_jis Px_jis_ss Pp_jis_ss M_ss K_ss Kbar_ss xpred_jis Ppredx_jis_ss Pxp_k_kmin_ss] = JIS_ss(mod.A,mod.B,mod.G,mod.J,sim.y,sim.x0,mod.Q,mod.R,mod.S,mod.P01,'trunc',false);

fid.L=10;
[x_smooth p_smooth Px_smooth_ss Pp_smooth_ss]=JIS_smooth(mod.A,mod.B,mod.G,mod.J,sim.y,mod.Q,mod.R,mod.S,sim.x0,mod.P01,fid.L,'convtol',1e-10);

plottime(sim.t,sim.p,p_jis_old,p_jis,p_smooth,'displayname',{'True' 'Filt old' 'Filt' 'Smooth'});
plotfreq(sim.t,sim.p,p_jis_old,p_jis,p_smooth,'displayname',{'True' 'Filt old' 'Filt' 'Smooth'});


% plotstat(diag(Pp_jis_ss_old),diag(Pp_jis_ss),plotopt);

plotopt=struct();
plotopt.gap=[0.1 0.1];
plotopt.marg_h=[0.1 0.25];
plotopt.ylabel='Error SD';
plotopt.markersize=10;
plotopt.color={[0 0 1] [1 0 0]   [0 0 1] [1 0 0]  };
plotopt.normalize=true;

% States
plot_std_emp_x=sdstatemp(sim.x,{xpred_jis , x_jis , x_smooth},[50 50])

plot_std_pred_x={};
plot_std_pred_x{1}=[diag(Ppredx_jis_ss).^0.5 ];
plot_std_pred_x{2}=[diag(Px_jis_ss).^0.5 ];
plot_std_pred_x{3}=[diag(Px_smooth_ss).^0.5 ];

plotstat(plot_std_emp_x{1},plot_std_pred_x{1},plotopt);
plotstat(plot_std_emp_x{2},plot_std_pred_x{2},plotopt);
plotstat(plot_std_emp_x{3},plot_std_pred_x{3},plotopt);

% Force
plot_std_emp_p=sdstatemp(sim.p,{[p_jis] , [p_smooth]},[50 50])

plot_std_pred_p{1}=[diag(Pp_jis_ss).^0.5];
plot_std_pred_p{2}=[diag(Pp_smooth_ss).^0.5];

plotstat(plot_std_emp_p{1},plot_std_pred_p{1},plotopt);

plotstat(plot_std_emp_p{2},plot_std_pred_p{2},plotopt);

tilefigs

return
%% Test lag

L_mat=[1 3 5 10 20 30 40 50];
% L_mat=[10:10:30];

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

tilefigs
%% KF

[x_k_k,x_k_kmin,x_k_N,P_k_k,P_k_kmin,P_k_N]=KF_RTSS(mod.A,mod.B,mod.G,mod.J,mod.Q,mod.R,mod.S,sim.y,sim.p);

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

%% Verify covariance between 

[x_jis p_jis Px_jis_ss Pp_jis_ss M_ss K_ss Kbar_ss x_pred,Px_k_kmin_jis_ss,Pxp_k_kmin_jis_ss] = JIS_ss(mod.A,mod.B,mod.G,mod.J,sim.y,sim.x0,mod.Q,mod.R,mod.S,mod.P01,'trunc',false);

P_emp=cov(([sim.x;sim.p]-[x_pred;p_jis]).');

P_th=[Px_k_kmin_jis_ss Pxp_k_kmin_jis_ss ; Pxp_k_kmin_jis_ss.' Pp_jis_ss];

dP=P_emp-P_th;

C_emp=cov2corr(P_emp)

C_th=cov2corr(P_th)

dC=C_th-C_emp

figure();
imagesc(dC,[-0.001 0.001]);
h = colorbar;
plotcorr(P_emp)
plotcorr(P_th)
tilefigs

close all

plotopt=struct();
plotopt.gap=[0.1 0.1];
plotopt.marg_h=[0.1 0.25];
plotopt.ylabel='Error SD';
% plotopt.displayname={'Pred' 'Filt' 'Smooth' };
plotopt.markersize=10;
plotopt.color={[0 0 1] [1 0 0] [0 0 0]  };
% plotopt.marker={ 'x' 'o' 'd'};
plotopt.normalize=true

plotstat(std([sim.x;sim.p]-[x_pred;p_jis],0,2),sqrt(diag(P_th)),plotopt);

tilefigs

e_jis=sim.y-mod.G*x_pred-mod.J*p_jis;

P_pred=[mod.G mod.J eye(size(mod.R))]*...
       [...
       Px_k_kmin_jis_ss Pxp_k_kmin_jis_ss zeros(size(Px_k_kmin_jis_ss,1),size(mod.R,1)) ; ...
       Pxp_k_kmin_jis_ss.' Pp_jis_ss  (-mod.R*M_ss.').' ; ...
       zeros(size(mod.R,1),size(Px_k_kmin_jis_ss,1)) (-mod.R*M_ss.') mod.R ...
       ]...
       *[mod.G mod.J eye(size(mod.R))].';


cov_e_jis=cov(e_jis.');

clc
close all
figure();
imagesc(dC,[-0.001 0.001]);
h = colorbar;
plotcorr(cov2corr(cov_e_jis))
plotcorr(cov2corr(P_pred))

tilefigs

close all

plotopt=struct();
plotopt.gap=[0.1 0.1];
plotopt.marg_h=[0.1 0.25];
% plotopt.normalize=true;
plotstat(diag(P_pred),diag(cov_e_jis),plotopt);


temp=cov([sim.p'-p_jis' sim.v.']);
temp2=temp(1,2:6)
check=(-mod.R*M_ss.').'

% 



% temp=cov([sim.p'-p_jis' sim.v.']);temp2=temp(6,1:5)




% P_pred_wrong=[mod.G mod.J eye(size(mod.R))]*...
%        [...
%        Px_k_kmin_jis_ss Pxp_k_kmin_jis_ss zeros(size(Px_k_kmin_jis_ss,1),size(mod.R,1)) ; ...
%        Pxp_k_kmin_jis_ss.' Pp_jis_ss  0*(-mod.R*M_ss.').' ; ...
%        zeros(size(mod.R,1),size(Px_k_kmin_jis_ss,1)) 0*(-mod.R*M_ss.') mod.R ...
%        ]...
%        *[mod.G mod.J eye(size(mod.R))].';

% e_st=sim.y-mod.G*x_pred-mod.J*p_jis;
% cov_e_st=cov(e_st.');
% Rk_st=mod.G*Px_k_kmin_jis_ss*mod.G.'+mod.R
% 
% Rk_st1=[mod.G mod.J]*...
%     [Px_k_kmin_jis_ss Pxp_k_kmin_jis_ss ; Pxp_k_kmin_jis_ss.' Pp_jis_ss]*...
%     [mod.G mod.J].'+mod.R




%%
clc
close all

Aa=[mod.A mod.B ; zeros(mod.nm*2,mod.np).' eye(mod.np)];

Ga=[mod.G mod.J]

Ra=mod.R;
Sa=[mod.S; zeros(mod.np,mod.ny)];
Qp=eye(mod.np)*1e2;
Qa=blkdiag(mod.Q,Qp);

[xa,~,~,P_k_k,~,~]=KF_RTSS(Aa,[],Ga,[],Qa,Ra,Sa,sim.y,[]);

xhat=xa(1:end-2,:);
phat=xa(end-1:end,:);

plottime(sim.t,sim.p,phat,'displayname',{'True' 'Filt'});


%%


% Check statistics (best for long time series, band limited white noise input)

% close all
% 
% plotopt=struct();
% plotopt.gap=[0.1 0.1];
% plotopt.marg_h=[0.1 0.25];
% plotopt.ylabel='Error SD';
% plotopt.displayname={'Pred' 'Filt' 'Smooth' };
% plotopt.markersize=10;
% plotopt.color={[0 0 1] [1 0 0] [0 0 0]  };
% plotopt.marker={ 'x' 'o' 'd'};


