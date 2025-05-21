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


%%


% omega_axis=linspace(0,20,1e4);

H=ssmod_tf(Fc,Lc,Hc,[],omega_axis);

S_w=sigma_w.^2/(2*pi)*ones(1,1,length(omega_axis));

S_par_3=mtimes3(H,S_w,H,'nnh');

S_par_3_onesided=S_par_3*2;

close all

plotpsd(omega_axis,S_par_3_onesided,w_welch,S_welch,...
    'xlabel','Frequency [rad/s]','ylabel','PSD',...
    'displayname',{'Test' 'Welch'},...
    'LineStyle',{'-' '--' '-'},...
    'log',false,'xlim',[0 5]);


return
%%



close all


% factor=1/dt^2;

% Exact covariances
Qc2=sigma_w^2
Qd2=Qc2/dt;
[Fd2,Ld2,Hd2,~]=ssmod_c2d(Fc,Lc,Hc,[],dt);


w2=mvnrnd(zeros(size(Fd,1),1),Qd2,length(t)).';

[s2,y2]=ssmod_forward(Fd2,Ld2,Hd2,[],[],[0],w2);

std_y2=std(y2)

plottime(t,y2,'ylabel','Simulated process');
xlim([0 1000])

% Empirical spectrum
[S_welch2,w_welch2]=estimateSpectrumWelch(y2,1/dt,'unit','rad','plot','no','Nwelch',100);


plotpsd(w_welch,S_welch,w_welch2,S_welch2,...
    'xlabel','Frequency [rad/s]','ylabel','PSD',...
    'LineStyle',{'-' '--' '-'},...
    'log',false,'xlim',[0 5]);


tilefigs



