%%


clc

[~,~,fid.Sp_prime]= generalSelection({},{},getLabel('U2',[60]+1e3),phi_label);

[~,mod.B_prime,~,mod.J_prime]=statespaceModel(mod.phi, mod.omega, mod.gamma, fid.Sa, fid.Sd, fid.Sp_prime, dt);

fid.p_prime=1/2*sin(2*fid.t)+1/3*sin(3*fid.t)+1/4*sin(4*fid.t);


[fid.x1 fid.y1]=statespaceForward(mod.A,mod.B,mod.G,mod.J,mod.F*0,fid.x0,fid.p);

[fid.x2 fid.y2]=statespaceForward(mod.A,mod.B_prime,mod.G,mod.J_prime,mod.F*0,fid.x0,fid.p_prime);

fid.x=fid.x1+fid.x2;

fid.y=fid.y1+fid.y2;

varR=var(fid.y,0,2)*0.0001^2;
fid.R=diag(varR);
fid.Q=eye(mod.nm*2)*1e-14;

fid.y_clean=fid.y;

rng(0); fid.y=fid.y_clean+randn(mod.nd,length(fid.t)).*varR.^0.5;

close all

plotTime(fid.t,fid.y1,fid.y2,fid.y);

plotTime(fid.t,fid.p,fid.p_prime);


%% Test: no compensation of deterministic forces

% close all
% 
% fid.L=30;
% [xsmooth psmooth P_x_ss P_p_ss ]=JIS_smooth(mod.A,mod.B,mod.G,mod.J,fid.y,fid.R,fid.Q,fid.S,fid.x0,fid.P01,fid.L);
% 
% [xfilt pfilt] = JIS_trunc(mod.A,mod.B,mod.G,mod.J,fid.y,fid.x0,fid.R,fid.Q,fid.S,fid.P01,'steadystate','yes');
% 
% plotTime(fid.t,fid.p,pfilt,psmooth,'legend',{'True' 'Filt' 'Smooth'});
% 
% plotFreq(fid.t,fid.p,pfilt,psmooth,'xlim',[0 50],'legend',{'True' 'Filt' 'Smooth'});
% 
% plotTime(fid.t,fid.x(1:6,:),xfilt(1:6,:),xsmooth(1:6,:),'legend',{'True' 'Filt' 'Smooth'});
% 
% plotFreq(fid.t,fid.x(1:6,:),xfilt(1:6,:),xsmooth(1:6,:),'xlim',[0 50],'legend',{'True' 'Filt' 'Smooth'});
% 
% tilefigs



%% Test: compensation of deterministic forces

clc
close all

[xfilt pfilt]=JIS_DI(mod.A,mod.B,mod.B_prime,mod.G,mod.J,mod.J_prime,fid.y,fid.p_prime,fid.x0,fid.R,fid.Q,fid.S,fid.P01);

plotTime(fid.t,fid.p,pfilt,'legend',{'True' 'Filt'});

plotFreq(fid.t,fid.p,pfilt,'xlim',[0 50],'legend',{'True' 'Filt' });

plotTime(fid.t,fid.x(1:6,:),xfilt(1:6,:),'legend',{'True' 'Filt'});

plotFreq(fid.t,fid.x(1:6,:),xfilt(1:6,:),'xlim',[0 50],'legend',{'True' 'Filt'});

tilefigs


