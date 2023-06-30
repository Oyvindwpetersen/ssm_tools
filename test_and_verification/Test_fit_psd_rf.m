%% Verification of fit of state-space models for
% for various parametric expressions
    
%% Single modal spectrum (exact rf)

clc
clear all
close all

omega=[0.01:0.01:20];

S_target(1,1,:)=400*(omega.^2)./(omega.^6-2*omega.^4+5*omega.^2+10);
order_n=2; order_d=6;
forcezero=true;

plotpsd(omega,S_target,'xlim',[0 20]);

S_target_twosided=S_target/2;
[n_opt,d_opt,alpha_opt]=fit_psd_rf(omega,S_target_twosided,order_n,order_d);

% plotpsd(omega,S_target_twosided,S_ssmod_twosided,'xlim',[0 10]);

%% Double modal spectrum (exponential)
    
clc
clear all
close all

omega=[0.01:0.01:10];

S_target(1,1,:)=20*exp(-0.5*(omega-2).^2/0.2^2)+10*exp(-0.5*(omega-4).^2/0.5^2);

forcezero=true; order_n=6; order_d=8;
forcezero=true; order_n=8; order_d=10;

% forcezero=false; order_n=2; order_d=4;
% forcezero=true; order_n=14; order_d=18;
[n_opt,d_opt,alpha_opt,S_opt,rn_opt,rd_opt]=fit_psd_rf(omega,S_target,order_n,order_d,'forcezero',forcezero);


%% Single modal spectrum with time domain sim

clc
clear all
close all

omega=[0.01:0.01:10];

omega_c=2.0
S_target(1,1,:)=omega.^-5.*exp(-5/4*omega_c.^4./omega.^4); S_target=S_target./max(S_target)*20; S_target(S_target<1e-12)=1e-12;

S_target_twosided=S_target/2;

forcezero=true; order_n=6; order_d=8;
% forcezero=true; order_n=4; order_d=6;
[n_opt,d_opt,alpha_opt,S_opt]=fit_psd_rf(omega,S_target_twosided,order_n,order_d,'forcezero',forcezero);

tilefigs([2 2],'l');

% State-space model
[Fc,Lc,Hc,sigma_w]=ssmod_psd_rf(n_opt,d_opt,alpha_opt);
Sw(1,1,1:length(omega))=sigma_w.^2/(2*pi);
Ht=ssmod_tf(Fc,Lc,Hc,zeros(1),omega);
S_ssmod_twosided=mtimes3(Ht,Sw,Ht,'nnh');

plotpsd(omega,S_target_twosided,S_opt,S_ssmod_twosided,'xlim',[0 10]);

% Simulate
dt=0.02;
Qc=Lc*sigma_w.^2*Lc.';
Qd=cov_c2d(Fc,Qc,dt);

[Fd,~,Hd,~]=ssmod_c2d(Fc,[],Hc,[],dt);
Bd=eye(size(Fd)); Jd=zeros(1,size(Hd,2));

t=[0:dt:100000];

w=mvnrnd(zeros(size(Qd,1),1),Qd,length(t)).';

x0=zeros(size(Fd,1),1);

[x,y]=ssmod_forward(Fd,Bd,Hd,Jd,[],x0,w);

plottime(t,y);

[S_welch,f_welch]=estimateSpectrumWelch(y,1/dt,'Nwelch',100); w_welch=f_welch*2*pi; S_welch=S_welch/(2*pi);

plotpsd(omega,S_target,omega,S_ssmod_twosided*2,w_welch,S_welch,'xlim',[0 10]);

