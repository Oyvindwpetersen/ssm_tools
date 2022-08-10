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

plotSpectrum(omega,S_exact,S_ss,'LineStyleSet',{'-' '--'},'xlim',[0 10],'legend',{'Exact original' 'State-space model'});
