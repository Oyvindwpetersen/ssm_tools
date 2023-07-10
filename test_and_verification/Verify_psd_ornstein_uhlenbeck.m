%% Verification of Ornstein Uhlenbeck process

clc
clear all
close all

dt=0.02

t=[0:dt:10000];
w_axis=[0:0.001:1/dt/2]*2*pi;

% Matern 0.5 
lambda=0.35;
Ac=-lambda; Bc=3; Gc=15; Jc=0;
[A B G J]=ssmod_c2d(Ac,Bc,Gc,Jc,dt);

% Bd becomes 1 or eye() when the system is discretized, see powerpoint
% Bc is taken into account in the covariance of the discrete noise 
B_sim=eye(1);

% White noise
sigma_wc_squared=8^2
sigma_wd_squared=Bc*sigma_wc_squared*Bc.'*dt;

S_eta_cont(1,1,1:length(w_axis))=sigma_wc_squared/(2*pi);
H_ss_cont=ssmod_tf(Ac,Bc,G,J,w_axis,[]);
S_y_ss_cont=mtimes3(H_ss_cont,S_eta_cont,H_ss_cont,'nnh');

H_ss_disc=ssmod_tf(A,B_sim,G,J,w_axis,dt);
S_eta_disc(1,1,1:length(w_axis))=sigma_wd_squared/(2*pi); % Not quite sure why this is, but it gives correct
S_y_ss_disc=mtimes3(H_ss_disc,S_eta_disc*dt,H_ss_disc,'nnh'); % Not quite sure why this is, but it gives correct

% Simulate in time
w_sim=randn(size(t))*sigma_wd_squared.^0.5;
[x_sim,y_sim]=ssmod_forward(A,B_sim,G,J,[],0,w_sim);
plottime(t,y_sim)

% From theory, see powerpoint
S_theory(1,1,:)=(Bc*Gc)^2*sigma_wc_squared/(2*pi)./(lambda.^2+w_axis.^2);

% Empirical from Welch
[S_welch,w_welch]=estimateSpectrumWelch(y_sim,1/dt,'Nwelch',100,'unit','rad','plot','no');

plotopt=struct()
plotopt.xlim=[0 100];
plotopt.linestyle={'-' '--' ':' '--'}

plotpsd(w_welch,S_welch,w_axis,S_y_ss_cont*2,w_axis,S_theory*2,plotopt);

tilefigs([],'l');
