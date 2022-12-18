%%

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
% mod.S=[ones(mod.nx/2,mod.ny)*0.1 ; ones(mod.nx/2,mod.ny)*-0.1].*diag(mod.Q).^0.5.*(diag(mod.R).').^0.5;

% Ctot=[mod.Q mod.S ; mod.S.' mod.R]
Ccorr=plotcorr([mod.Q mod.S ; mod.S.' mod.R]);

sim.nt=1e5;
sim.t=[1:sim.nt]*mod.dt;

%% Generate noise

close all
% sim.p=BandLimitedWhiteNoise(0.1,3,sim.t,mod.np)*1;
sim.p=BellShapedNoise(0.1,1,sim.t,mod.np)*1;

[sim.w,sim.v]=cov_noisegen(mod.Q,mod.R,mod.S,sim.t);

sim.x0=zeros(mod.nx,1);
[sim.x,sim.y]=ssmod_forward_stoch(mod.A,mod.B,mod.G,mod.J,[],sim.x0,sim.p,sim.w,sim.v);

plotTime(sim.t,sim.x);
plotTime(sim.t,sim.y);

tilefigs


%% Estimate

close all

mod.P01=eye(size(mod.Q));

[x_jis p_jis Px_jis_ss Pp_jis_ss ] = JIS_trunc_ss(mod.A,mod.B,mod.G,mod.J,sim.y,sim.x0,mod.R,mod.Q,mod.S,mod.P01,'trunc','no');

fid.L=60;
[x_smooth p_smooth Px_smooth_ss Pp_smooth_ss ]=JIS_smooth(mod.A,mod.B,mod.G,mod.J,sim.y,mod.R,mod.Q,mod.S,sim.x0,mod.P01,fid.L);

plotTime(fid.t,fid.p,p_jis,p_smooth,'legend',{'True' 'Filt' 'Smooth'});

plotFreq(fid.t,fid.p,p_jis,p_smooth,'xlim',[0 50],'legend',{'True' 'Filt' 'Smooth'});

plotTime(fid.t,fid.x(1:6,:),x_jis(1:6,:),x_smooth(1:6,:));

plotFreq(fid.t,fid.x(1:6,:),x_jis(1:6,:),x_smooth(1:6,:),'xlim',[0 50]);

tilefigs

%% Check statistics (best for long time series, band limited white noise input)

err_x_jis=std(fid.x-x_jis,0,2) %./ std(fid.x,0,2);
err_x_smooth=std(fid.x-x_smooth,0,2) %./ std(fid.x,0,2);

err_p_jis=std(fid.p-p_jis,0,2) %./ std(fid.p,0,2);
err_p_smooth=std(fid.p-p_smooth,0,2) %./ std(fid.p,0,2);

close all

% figure(); hold on; grid on; ylog;
% plot(err_x_jis,'ob','DisplayName','Error JIS');
% plot(err_x_smooth,'xr','DisplayName','Error Smooth');
% plot(diag(Px_jis_ss).^0.5,'*k','DisplayName','Error JIS (from filter)');
% plot(diag(Px_smooth_ss).^0.5,'dk','DisplayName','Error Smooth (from filter)');
% legend show
% 
% figure(); hold on; grid on; ylog;
% plot(err_p_jis,'ob','DisplayName','Error JIS');
% plot(err_p_smooth,'xr','DisplayName','Error Smooth');
% plot(diag(Pp_jis_ss).^0.5,'*k','DisplayName','Error JIS (from filter)');
% plot(diag(Pp_smooth_ss).^0.5,'dk','DisplayName','Error Smooth (from filter)');
% legend show


figure(); hold on; grid on; ylog;
title('Ratio variance state'); legend show;
plot(err_x_jis./diag(Px_jis_ss).^0.5,'ob','DisplayName',' JIS');
plot(err_x_smooth./diag(Px_smooth_ss).^0.5,'*r','DisplayName',' Smooth');


figure(); hold on; grid on; ylog;
title('Ratio variance'); legend show;
plot(err_p_jis./diag(Pp_jis_ss).^0.5,'ob','DisplayName',' JIS');
plot(err_p_smooth./diag(Pp_smooth_ss).^0.5,'*r','DisplayName',' Smooth');

% 
% figure(); hold on; grid on; ylog;
% plot(err_p_jis,'ob','DisplayName','Error JIS');
% plot(err_p_smooth,'xr','DisplayName','Error Smooth');
% plot(diag(Pp_jis_ss).^0.5,'*k','DisplayName','Error JIS (from filter)');
% plot(diag(Pp_smooth_ss).^0.5,'dk','DisplayName','Error Smooth (from filter)');
% legend show
% 
% 

tilefigs


%% Test lag

% L_mat=[1 3 5 10 20 30 40 ];
L_mat=[10:10:50];

for k=1:length(L_mat)
    [~,~,P_x_ss,P_p_ss ]=JIS_smooth(mod.A,mod.B,mod.G,mod.J,sim.y,mod.R,mod.Q,mod.S,sim.x0,mod.P01,L_mat(k));
    Px_all(k,:)=diag(P_x_ss);
    Pp_all(k,:)=diag(P_p_ss);
end

close all
figure(); hold on; grid on;
plot(L_mat,Px_all); ylog;
    
figure(); hold on; grid on;
plot(L_mat,Pp_all); ylog;
