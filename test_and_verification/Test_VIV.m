%%

clc
clear all
close all

%% Create system

mod=importbeamquick(6);

mod.dt=0.02;

mod.a_cell={'10_U'  '50_U'  '80_U' }
mod.d_cell={'10_U'  '50_U'  '80_U' }
mod.p_cell={}

[mod.Sd,mod.Sa,mod.Sp]=DofSelection(mod.d_cell,mod.a_cell,mod.p_cell,mod.doflabel);

[mod.A mod.B mod.G mod.J mod.Ac mod.Bc mod.Gc mod.Jc]=ssmod_modal(mod.phi,mod.Omega,mod.Gamma,mod.Sa,mod.Sd,mod.Sp,mod.dt,'force','modal');

mod.ny=size(mod.Sa,1); mod.nx=size(mod.A,1); mod.np=size(mod.Sp,2);

%% Load

sim.nt=1e5+1;
sim.t=[1:sim.nt]*mod.dt;

sigma_buff=80

sim.p_buff=randn(6,sim.nt)*sigma_buff;

% sim.p_viv=cos(2.7*sim.t)*1000;

rng(0);
sim.p_viv=bandlimitedwn(0.44,0.46,sim.t,1)*1000;
sim.p_viv=sim.p_viv+randn(size(sim.p_viv));

plottime(sim.t,sim.p_viv);

%% Simulate

close all

Bc_viv=mod.B(:,2);
Bc_buff=mod.B(:,[1 2 3 4 5 6]);

Jc_viv=mod.J(:,2);
Jc_buff=mod.J(:,[1 2 3 4 5 6]);

[mod.Ad mod.Bd_viv mod.Gd mod.Jd_viv]=ssmod_c2d(mod.Ac,Bc_viv,mod.Gc,Jc_viv,mod.dt);
[mod.Ad mod.Bd_buff mod.Gd mod.Jd_buff]=ssmod_c2d(mod.Ac,Bc_buff,mod.Gc,Jc_buff,mod.dt);

R0=1e-8*eye(mod.ny);
v_noise=mvnrnd(zeros(size(R0,1),1),R0,sim.nt).';

w=mod.Bd_buff*sim.p_buff;
v=mod.Jd_buff*sim.p_buff+v_noise;

[sim.x,sim.y]=ssmod_forward_stoch(mod.Ad,mod.Bd_viv,mod.Gd,mod.Jd_viv,[],0,sim.p_viv,w,v);


plottime(sim.t,sim.x);
plotfreq(sim.t,sim.x,'xlim',[0 10]);

plottime(sim.t,sim.y);
plotfreq(sim.t,sim.y,'xlim',[0 10]);

tilefigs

%% Estimate

Q=mod.Bd_buff*sigma_buff^2*eye(6)*mod.Bd_buff.';
R=mod.Jd_buff*sigma_buff^2*eye(6)*mod.Jd_buff.'+R0;
S=mod.Bd_buff*sigma_buff^2*eye(6)*mod.Jd_buff.';

[x_jis p_jis Px_jis_ss Pp_jis_ss] = JIS_ss(mod.Ad,mod.Bd_viv,mod.Gd,mod.Jd_viv,sim.y,0,Q,R,S,[]);

[x_smooth p_smooth Px_smooth_ss Pp_smooth_ss] = JIS_smooth(mod.Ad,mod.Bd_viv,mod.Gd,mod.Jd_viv,sim.y,Q,R,S,0,[],20);

[A1 B1 G1 J1]=ssmod_modal(mod.phi(:,2),mod.Omega(2,2),mod.Gamma(2,2),mod.Sa,mod.Sd,mod.Sp,mod.dt,'force','modal');

[x_jis_1 p_jis_1 Px_jis_ss Pp_jis_ss] = JIS_ss(A1,B1,G1,J1,sim.y,0,zeros(2),R0,zeros(2,6),[]);

% [x_smooth p_smooth Px_jis_ss Pp_jis_ss] = JIS_smooth(mod.Ad,mod.Bd_viv,mod.Gd,mod.Jd_viv,sim.y,Q,R,S,0,[],20);

x_jis_1_resize=[zeros(1,sim.nt) ; x_jis_1(1,:) ; zeros(4,sim.nt) ; ...
    zeros(1,sim.nt) ; x_jis_1(2,:) ; zeros(4,sim.nt) ; ...
    ];


%% Plot

close all

plotopt=struct();
plotopt.linestyle={'-' '--' ':'};

plottime(sim.t,sim.x,x_jis,x_jis_1_resize,plotopt);
plotfreq(sim.t,sim.x,x_jis,x_jis_1_resize,plotopt);

plottime(sim.t,sim.p_viv,p_jis,p_jis_1,plotopt);


plotfreq(sim.t,sim.p_viv,p_jis,p_jis_1,plotopt,'xlim',[0 50]);

tilefigs


