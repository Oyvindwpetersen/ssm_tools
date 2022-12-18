%%

% fid.R=diag([1e-3 1e-3 1e-4 1e-4].^2]);
% % fid.noise=mvnrnd(zeros(size(fid.y_clean,1),1),fid.R,length(fid.t)).';
% fid.y=fid.y_clean+fid.noise;
% 
% fid.Q=eye(mod.nm*2)*1e-12;
% fid.S=zeros(mod.nm*2,mod.ny);
% 
% fid.P01=fid.Q;
% fid.x0=zeros(mod.nm*2,1);

%%
clc
close all

lambda=0.05;
sigma_p1=1e-3;
[Fc1,Lc1,Hc1,sigma_w1]=ssmod_matern(lambda,0,sigma_p1)
[Fd1,Ld1,Hd1]=ssmod_c2d(Fc1,Lc1,Hc1,[],dt);

lambda=0.05;
sigma_p2=1e-5;
[Fc2,Lc2,Hc2,sigma_w2]=ssmod_matern(lambda,0,sigma_p2)
[Fd2,Ld2,Hd2]=ssmod_c2d(Fc2,Lc2,Hc2,[],dt);

Fc=diag([Fc1;Fc1;Fc2;Fc2]);
Lc=diag([Lc1;Lc1;Lc2;Lc2]);

Fd=diag([Fd1;Fd1;Fd2;Fd2]);
Ld=diag([Ld1;Ld1;Ld2;Ld2]);
Hd=diag([Hd1;Hd1;Hd2;Hd2]);

Sigmac=diag([sigma_w1;sigma_w1;sigma_w2;sigma_w2].^2);
Sigmad=cov_c2d(Fc,Lc*Sigmac*Lc.',dt);

vbar=mvnrnd(zeros(4,1),Sigmad,length(fid.t)).';
[~,v_col_sim]=ssmod_forward(Fd,eye(4),Hd,zeros(4,4),[],zeros(4,1),vbar);

[R_emp tau_emp]=xcorrSignal(v_col_sim,dt,10000);
tau_theory=[0:0.1:50];
R_theory=matern_cf(tau_theory,[lambda lambda lambda lambda],[sigma_w1 sigma_w1 sigma_w2 sigma_w2],[0 0 0 0],'3d');
plotAutocorr(tau_emp,R_emp,tau_theory,R_theory,'xlim',[0 100]);

% return
Zd=diag([1e-3 1e-3 1e-5 1e-5].^2);

% Zd=Zc*dt;

v_wn_sim=mvnrnd(zeros(4,1),Zd,length(fid.t)).';

plotTime(fid.t,v_wn_sim,v_col_sim);

[S_welch,f_welch]=estimateSpectrumWelch(v_wn_sim+v_col_sim,1/dt);

plotSpectrum(f_welch,S_welch,'xlim',[0 10]);

%%
fid.y=fid.y_clean+v_wn_sim+v_col_sim;

% 
% 

fid.R=Zd+Sigmad;
fid.Q=eye(size(mod.A))*1e-12;
fid.S=zeros(mod.nm*2,mod.ny);

fid.P01=fid.Q;
fid.x0=zeros(mod.nm*2,1);

close all

[x_jis p_jis Px_jis_ss Pp_jis_ss ] = JIS_trunc_ss(mod.A,mod.B,mod.G,mod.J,fid.y,fid.x0,fid.R,fid.Q,fid.S,fid.P01,'trunc','yes');

plotTime(fid.t,fid.p,p_jis);


return

%%

mod2=struct();
fid2=struct();

mod2.A=blkdiag(mod.A,Fd);
mod2.B=[mod.B ; zeros(4,1)];

mod2.G=[mod.G Hd];
mod2.J=mod.J;

fid2.Q=blkdiag(fid.Q,Sigmad);

fid2.R=Zd;
fid2.S=zeros(mod.nm*2+4,mod.ny);

fid2.P01=fid2.Q;
fid2.x0=zeros(mod.nm*2+4,1);


[x_jis2 p_jis2 Px_jis_ss Pp_jis_ss ] = JIS_trunc_ss(mod2.A,mod2.B,mod2.G,mod2.J,fid.y,fid2.x0,fid2.R,fid2.Q,fid2.S,fid2.P01,'trunc','yes');


close all


plotTime(fid.t,fid.p,p_jis,p_jis2);
%%

close all

plotFreq(fid.t,fid.p-p_jis,fid.p-p_jis2);


% plotTime(fid.t,v_col_sim);


% 
% w2=mvnrnd(zeros(4,1),Q_sim,length(t)).';
% v=mvnrnd(zeros(4,1),R_sim,length(t)).';


