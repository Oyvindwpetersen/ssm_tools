%%

clc
clear all
close all

%% Model

load('abaqus\ssbeam_twospan.mat')

mod.includeModes=[1:6];

mod.omega=diag(freq(mod.includeModes)*2*pi);
mod.xi=eye(size(mod.omega))*0.002;
mod.gamma=2*mod.xi.*mod.omega;
mod.phi=phi(:,mod.includeModes);


%% Plot modes
close all

figure(); hold on;

plot(mod.phi(2:6:end,1:6));
% plot(mod.phi(2:6:end,1));

linevertical(gca,[20 40],'-',2,'b')
linevertical(gca,[10],'-',2,'k')

%% Forces

mod.a_cell=genlabels('U2',[20 40]+1e3); 
mod.p_cell=genlabels('U2',[12 ]+1e3);
mod.d_cell=genlabels('U2',[20 40]+1e3);


mod.a_cell=genlabels('U2',[20 40 ]+1e3); 
mod.p_cell=genlabels('U2',[20 ]+1e3);
mod.d_cell=genlabels('U2',[20 40]+1e3);

% mod.a_cell=genlabels('U2',[10 20 40 50 60 70 ]+1e3); 
% mod.p_cell=genlabels('U2',[10 20 40 50 60 70 ]+1e3);
% mod.d_cell=genlabels('U2',[10 20 40 50 60 70 ]+1e3);

[mod.Sd, mod.Sa, mod.Sp]= DofSelection(mod.d_cell,mod.a_cell,mod.p_cell,phi_label);

dt=100^-1;
Fs=dt^-1;
[mod.A mod.B mod.G mod.J mod.Ac mod.Bc mod.F]=ssmod_modal(mod.phi, mod.omega, mod.gamma, mod.Sa, mod.Sd, mod.Sp,dt,'force','disc');

mod.nm=size(mod.phi,2)
mod.ny=length(mod.a_cell)+length(mod.d_cell);
mod.np=length(mod.p_cell);
% mod.np=mod.nm;

% Check requirements
tz=ssmod_tzero(mod.A,mod.B,mod.G,mod.J);


%%

T=1000
fid.t=[0:dt:T];

fid.x0=zeros(mod.nm*2,1);

fid.p=[];
% fid.p(1,:)=sin( cos(fid.t*1.5).*(fid.t+20)*0.2);
fid.p(1,:)=sin( cos(fid.t.^(1*0.5)).*fid.t*1/10);


fid.p(1,:)=chirp(fid.t-T/2,0.1,T/2,8,'quadratic',[],'convex');

pspectrum(fid.p,fid.t,'spectrogram','TimeResolution',1, ...
    'OverlapPercent',90,'Leakage',0.85)
% fid.p=BandLimitedWhiteNoise(0.4,8,fid.t,1)

% fid.p(1,:)=sin( 2*pi*linspace(0,10,length(fid.t)).*fid.t);
% fid.p(1,:)=sin( 2*pi*(exp(linspace(-2,1,length(fid.t))/2)).*fid.t);

% fid.p_d=zeros(6,length(fid.t));
% for k=1:6
% fid.p_d(k,:)=sin( sin(fid.t.^(k*0.2)).*fid.t*0.5);
% end
% fid.p=mod.phi.'*mod.Sp*fid.p_d;

[fid.x fid.y_clean]=ssmod_forward(mod.A,mod.B,mod.G,mod.J,mod.F*0,fid.x0,fid.p);

close all

plottime(fid.t,fid.x);
plottime(fid.t,fid.y_clean);
plotfreq(fid.t,fid.y_clean,'xlim',[0 10]);

plottime(fid.t,fid.p);
plotfreq(fid.t,fid.p,'xlim',[0 10]);

tilefigs;

%%

fid.R=diag(std(fid.y_clean,0,2).^2)*0.01^2;
fid.noise=mvnrnd(zeros(size(fid.y_clean,1),1),fid.R,length(fid.t)).';
fid.y=fid.y_clean+fid.noise;

fid.Q=eye(mod.nm*2)*1e-12;
fid.S=zeros(mod.nm*2,mod.ny);

fid.P01=fid.Q;
fid.x0=zeros(mod.nm*2,1);


%%