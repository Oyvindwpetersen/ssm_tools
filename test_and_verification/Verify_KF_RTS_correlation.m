%%

clc
close all

% Estimate by adjusting both KF and RTS
[x_k_k,x_k_kmin,x_k_N,P_k_k,P_k_kmin,P_k_N]=...
    KF_RTSS(mod.A,mod.B,mod.G,mod.J,mod.Q,mod.R,mod.S,sim.y,sim.p,'method','standard');


% Estimate by adjusting both KF and RTS
[x_k_k0,x_k_kmin0,x_k_N0,P_k_k0,P_k_kmin0,P_k_N0]=...
    KF_RTSS(mod.A,mod.B,mod.G,mod.J,mod.Q,mod.R,mod.S,sim.y,sim.p,'method','decorr');

% % % Run KF unadjusted
% [x_k_k,x_k_kmin,P_k_k,P_k_kmin]=...
%     KF(mod.A,mod.B,mod.G,mod.J,mod.Q,mod.R,mod.S,sim.y,sim.p,[],[]);
%
% % Run adjusted RTS
% [x_k_N,P_k_N]=...
%     RTSS(mod.A-mod.S/mod.R*mod.G*0.0000,x_k_k,x_k_kmin,P_k_k,P_k_kmin);

% [x_k_N,P_k_N]=...
%     RTSS(eye(size(mod.A))*0.99,x_k_k,x_k_kmin,P_k_k,P_k_kmin);


% Check that they are equal
% plottime(sim.t,x_k_kmin,x_k_kmin0);
plottime(sim.t,x_k_k,x_k_kmin0);
plottime(sim.t,x_k_N,x_k_N0);
plotfreq(sim.t,x_k_N,x_k_N0);

plottime(sim.t,x_k_N-x_k_N0);

ratio1=P_k_kmin./P_k_kmin0
ratio2=P_k_k./P_k_k0
ratio3=P_k_N./P_k_N0

tilefigs

plottime(sim.t,sim.x,x_k_k,x_k_N);


%%

close all

plotopt=struct();
plotopt.gap=[0.1 0.1];
plotopt.marg_h=[0.1 0.25];
plotopt.ylabel='Error SD';
plotopt.displayname={'Emp' 'Pred' };
plotopt.markersize=10;
plotopt.color={[0 0 1] [1 0 0] [0 0 0]   [0 0 1] [1 0 0]  [0 0 0]};
plotopt.marker={ 'x' 'o' };
plotopt.spread=false;
plotopt.normalize=true;

% States
plot_std_emp_x=sdstatemp(sim.x,{x_k_kmin , x_k_k , x_k_N},[50 50])

plot_std_pred_x={};
plot_std_pred_x{1}=diag(P_k_kmin).^0.5;
plot_std_pred_x{2}=diag(P_k_k).^0.5;
plot_std_pred_x{3}=diag(P_k_N).^0.5;

plotstat(plot_std_emp_x{1},plot_std_pred_x{1},plotopt);
plotstat(plot_std_emp_x{2},plot_std_pred_x{2},plotopt);
plotstat(plot_std_emp_x{3},plot_std_pred_x{3},plotopt);

tilefigs
