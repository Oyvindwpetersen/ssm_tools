%%


clc
clear all
close all


phi1=0.8

sigma_eps=2

dt=0.1;
t_step=[1:1000000];
t=t_step*dt;

epsilon=randn(1,length(t_step))*sigma_eps;

x=0
for k=1:(length(t)-1)
    x(k+1)=phi1*x(k)+epsilon(k);
end

sigma_x=std(x)


plotTime(t,x);

[S_welch,w_welch]=estimateSpectrumWelch(x,1/dt,'Nwelch',1000,'unit','rad');


w_axis=[0:0.01:5]*2*pi;
S_AR(1,1,:)=AR_PSD(w_axis,phi1,sigma_eps,1,dt);


% plotSpectrum(w_welch,S_welch);

plotSpectrum(w_welch,S_welch,w_axis,S_AR*2,'xlim',[0 30]);


