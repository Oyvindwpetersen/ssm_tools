%%

clc
clear all
close all

%%
close all

rng(1);

dt=0.01;

sigma_p=200
lambda=0.1
p_order=0
[~,~,~,sigma_w]=ssmod_matern(lambda,p_order,sigma_p)

% [Fc,Lc,Hc]=ssmod_matern(lambda,p_order)

omega_axis=[0:1e-3:10];
tau_axis=[0:0.1:100];

[Fc,Lc,Hc,Hct,S_par_1,R_par_1]=ssmod_maternglobal(lambda,sigma_w,p_order,omega_axis,tau_axis,true);
[~,~,~,~,S_par_2,R_par_2]=ssmod_maternglobal(lambda,sigma_w,p_order,omega_axis,tau_axis,false);


plotpsd(omega_axis,S_par_1,omega_axis,S_par_2,...
    'xlabel','Frequency [rad/s]','ylabel','PSD',...
    'displayname',{'Fast calculation' 'Slow calculation'},...
    'LineStyle',{'-' '--'},...
    'log',false,'xlim',[0 5]);

plottime(tau_axis,R_par_1,tau_axis,R_par_2,...
    'xlabel','Time lag [s]','ylabel','Autocorr',...
    'displayname',{'Fast calculation' 'Slow calculation'},...
    'LineStyle',{'-' '--'},...
    'log',false,'xlim',[0 50]);


tilefigs


%% Simulate time series, compare spectra

close all

dt=0.25;
t=[0:dt:1e5];

% Exact covariances
Qc=Lc*sigma_w^2*Lc.'
Qd=cov_c2d(Fc,Qc,dt);
[Fd,Ld_notused,Hd,~]=ssmod_c2d(Fc,eye(size(Fc)),Hc,[],dt);
Ld=eye(size(Ld_notused));
w=mvnrnd(zeros(size(Fd,1),1),Qd,length(t)).';
Dd=zeros(1,size(Hd,2));

[s,y]=ssmod_forward(Fd,Ld,Hd,Dd,[],zeros(size(Fd,1),1),w);

std_y=std(y)

plottime(t,y,'ylabel','Simulated process');
xlim([0 1000])

% Empirical spectrum
[S_welch,w_welch]=estimateSpectrumWelch(y,1/dt,'unit','rad','plot','no','Nwelch',100);

% Onesided analytical spectrum
S_par_1_onesided=2*S_par_1;
S_par_2_onesided=2*S_par_2;

plotpsd(omega_axis,S_par_1_onesided,omega_axis,S_par_2_onesided,w_welch,S_welch,...
    'xlabel','Frequency [rad/s]','ylabel','PSD',...
    'displayname',{'Fast calculation' 'Slow calculation' 'Welch'},...
    'LineStyle',{'-' '--' '-'},...
    'log',false,'xlim',[0 5]);


tilefigs





