%% Verification of approximation of SE kernel 

clc
clear all
close all

rng(1);

% Hyperparameters
sigma_p=3.5
L=1.2
ns=4;

% Approximate state space model
omega_axis=[0:0.01:10];
[Fc,Lc,Hc,sigma_w,S_trunc]=ssmod_squaredexp(L,sigma_p,ns,omega_axis); S_trunc_onesided=S_trunc*2;

%% Simulate time series, compare spectra

close all

dt=0.025;
t=[0:dt:1e5];

% Exact covariances
Qc=Lc*sigma_w^2*Lc.'
Qd=cov_c2d(Fc,Qc,dt);
[Fd,~,Hd,~]=ssmod_c2d(Fc,[],Hc,[],dt);
Ld_sim=eye(size(Qd));

% Approximated covariances
% Qc=sigma_w^2;
% Qd=Qc*dt;
% [Fd,Ld,Hd,~]=ssmod_c2d(Fc,Lc,Hc,[],dt);

w=mvnrnd(zeros(size(Fd,1),1),Qd,length(t)).';
Dd=zeros(1,size(Hd,2));

[s,y]=ssmod_forward(Fd,Ld_sim,Hd,Dd,[],zeros(size(Fd,1),1),w);

std_y=std(y)

plottime(t,y,'ylabel','Simulated process');
xlim([0 1000])

[S_welch,w_welch]=estimateSpectrumWelch(y,1/dt,'unit','rad','plot','no','Nwelch',100);

% Exact PSD for SE
S_SE_exact(1,1,:)=exp(-0.5*L^2*omega_axis.^2) * sigma_p^2*L/sqrt(2*pi);
S_SE_exact_onesided=S_SE_exact*2;
sigma_p_test=trapz(omega_axis,S_SE_exact_onesided).^0.5

% PSD for state space model
Ht=ssmod_tf(Fc,Lc,Hc,zeros(1),omega_axis,[],'type','io');
S_w=repmat(sigma_w.^2/(2*pi),1,1,length(omega_axis));
S_ss_onesided=mtimes3(Ht,S_w,Ht,'nnh')*2;

plotpsd(omega_axis,S_SE_exact_onesided,w_welch,S_welch,omega_axis,S_ss_onesided,omega_axis,S_trunc_onesided,...
    'xlabel','Frequency [rad/s]','ylabel','PSD',...
    'displayname',{'SE onesided (exact)' 'Empirical (welch)' 'State-space onesided' 'Poly Taylor onesided'},...
    'LineStyle',{'-' '-' '--' ':'},...
    'log',false,'xlim',[0 10]);

tilefigs

%% Plot CF

tau_axis=[0:dt:100];

clear R_SE_exact
R_SE_exact(1,1,:)=sigma_p^2*exp(-0.5*tau_axis.^2/L^2);

[R_emp tau_emp]=xcorrSignal(y,dt,10000);

% close all

plottime(tau_axis,R_SE_exact,tau_emp,R_emp,...
    'xlabel','Time lag [s]','ylabel','ACF',...
    'displayname',{'R SE (exact)'  'Empirical (xcorr)'},...
    'LineStyle',{'-' '-' '--' ':'},'xlim',[0 50]);

tilefigs
