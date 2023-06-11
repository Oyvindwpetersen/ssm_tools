%%

clc
clear all
close all

omega=[0:0.01:20];

a=[1 -3 5 2];
n=[0 1 0];
S_exact(1,1,:)=omega.^2./(1-3*omega.^2+5*omega.^4+2*omega.^6);

[Fc,Lc,Hc,sigma_w,alpha]=ssmod_psd_rf(n,a);

H=ssmod_tf(Fc,Lc,Hc,zeros(1,1),omega);

Sw=ones(1,1,length(omega))*sigma_w.^2/(2*pi);

S_ss=mtimes3(H,Sw,H,'nnh');

close all

%% Time domain sim

dt=0.02;
[Fd,~,Hd]=ssmod_c2d(Fc,[],Hc,[],dt);
Ld=eye(size(Fd));

Qc=Lc*sigma_w.^2*Lc.';
Qd=cov_c2d(Fc,Qc,dt);

t=[0:dt:10000];
s0=zeros(size(Fd,1),1);

w_sim=mvnrnd(zeros(size(Fd,1),1),Qd,length(t)).';

[s,p]=ssmod_forward(Fd,Ld,Hd,zeros(size(Fd,1),1).',[],s0,w_sim);

[S_welch_temp,f_welch]=estimateSpectrumWelch(p,1/dt,'Nwelch',10);

w_welch=f_welch*2*pi;
S_welch=S_welch_temp/(2*pi)*0.5;

close all

plotTime(t,p);

%%
close all

plotpsd(omega,S_exact,...
    omega,S_ss,...
    w_welch.',S_welch,...
    'LineStyleSet',{'-' '--' '-'},'xlim',[0 10],'legend',{'Exact original' 'State-space model' 'Simulation'});



