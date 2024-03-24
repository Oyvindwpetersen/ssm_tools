%% Test PSD of WN

clc
clear all
close all

dt=0.02

t=[0:dt:100000];
omega_axis=[0:0.001:1/dt/2]*2*pi;

sigma_wc_squared=8^2

sigma_wd_squared=sigma_wc_squared*dt;

w_sim=randn(size(t))*sigma_wd_squared.^0.5;

[S_welch,w_welch]=estimateSpectrumWelch(w_sim,1/dt,'Nwelch',1000,'unit','rad');

S_test(1,1,1:length(omega_axis))=sigma_wc_squared/(2*pi);

plotpsd(w_welch,S_welch,omega_axis,S_test*2*dt^2,'xlim',[0 500]); 

% Still not quite sure why dt^2 is needed here, but ok
% Probably related to energy within Nyquist range
