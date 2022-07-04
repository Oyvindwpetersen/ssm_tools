%%

clc
clear all
close all

%%

rng(1);

dt=0.01;

lambda=0.1
sigma_w=2
p_order=2

% [Fc,Lc,Hc]=ssmod_matern(lambda,p_order)

w_axis=[0:1e-3:10];
tau_axis=[0:0.1:100];

[Fc,Lc,Hc,Ht,S_par_1,R_par_1]=ssmod_maternglobal(lambda,sigma_w,p_order,w_axis,tau_axis,true);
[~,~,~,~,S_par_2,R_par_2]=ssmod_maternglobal(lambda,sigma_w,p_order,w_axis,tau_axis,false);


plotSpectrum(w_axis,S_par_1,w_axis,S_par_2,...
    'xlabel','Frequency [rad/s]','ylabel','PSD',...
    'legend',{'Fast calculation' 'Slow calculation'},...
    'LineStyle',{'-' '--'},...
    'log','no','xlim',[0 5]);


plotAutocorr(tau_axis,R_par_1,tau_axis,R_par_2,...
    'xlabel','Time lag [s]','ylabel','Autocorr',...
    'legend',{'Fast calculation' 'Slow calculation'},...
    'LineStyle',{'-' '--'},...
    'log','no','xlim',[0 50]);



%% Simulate time series, compare spectra

close all

dt=0.25;
t=[0:dt:1e5];

% Exact covariances
Qc=Lc*sigma_w^2*Lc.'
Qd=cov_c2d(Fc,Qc,dt);
[Fd,Ld,Hd,~]=ssmod_c2d(Fc,eye(size(Fc)),Hc,[],dt);
Ld=eye(size(Ld));
w=mvnrnd(zeros(size(Fd,1),1),Qd,length(t)).';
Dd=zeros(1,size(Hd,2));

[s,y]=ssmod_forward(Fd,Ld,Hd,Dd,[],zeros(size(Fd,1),1),w);

std_y=std(y)

plotTime(t,y,'ylabel','Simulated process');
xlim([0 1000])

[S_welch,w_welch]=estimateSpectrumWelch(y,1/dt,'unit','rad','plot','no','Nwelch',1000);

S_par_1_onesided=2*S_par_1;
S_par_2_onesided=2*S_par_2;

plotSpectrum(w_axis,S_par_1_onesided,w_axis,S_par_2_onesided,w_welch,S_welch,...
    'xlabel','Frequency [rad/s]','ylabel','PSD',...
    'legend',{'Fast calculation' 'Slow calculation' 'Welch'},...
    'LineStyle',{'-' '--' '-'},...
    'log','no','xlim',[0 5]);


tilefigs






