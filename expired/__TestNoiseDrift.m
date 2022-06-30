%%


fid.t=[0:dt:50];

fid.x0=zeros(mod.nm*2,1);

fid.p=[];

fid.p(1,:)=sin( sin(fid.t*1.5).*fid.t);
% fid.p(2,:)=cos( sin(fid.t*1.5).*(fid.t-20));

close all
plotTime(fid.t,fid.p);
plotFreq(fid.t,fid.p,'xlim',[0 50]);

tilefigs


[fid.x fid.y]=statespaceForward(mod.A,mod.B,mod.G,mod.J,mod.F*0,fid.x0,fid.p);
fid.y_clean=fid.y;


plotTime(fid.t,fid.y);
plotFreq(fid.t,fid.y,'xlim',[0 50]);

tilefigs

%%



noise.lambda=0.00001;
noise.F=exp(-dt*noise.lambda)*eye(mod.ny);

noise.L=eye(mod.ny);
noise.H=eye(mod.ny);

sigma_noise=std(fid.y_clean,0,2)*0.01;
noise_input=sigma_noise.*randn(size(noise.F,1),length(fid.t));

[s,u]=statespaceForward(noise.F,noise.L,noise.H,zeros(size(noise.F)),[],zeros(size(noise.F,1),1),noise_input);


close all
plotTime(fid.t,u);
plotFreq(fid.t,u,'xlim',[0 50]);

fid.y=fid.y_clean+u;

close all
plotTime(fid.t,fid.y_clean,fid.y);
plotFreq(fid.t,fid.y_clean,fid.y,'xlim',[0 50]);

tilefigs


%% Model with white noise

[fid.Q]=eye(size(mod.A))*1e-10;
[fid.R]=diag(std(u,0,2)).^2;
fid.S=zeros(mod.nm*2,mod.ny);

fid.P01=fid.Q;
fid.x0=zeros(mod.nm*2,1);

fid.p0=zeros(mod.np,1);
fid.Pp0=10*eye(mod.np);
fid.Qp=1e3*eye(mod.np);

[xfilt pfilt] = JIS_trunc(mod.A,mod.B,mod.G,mod.J,fid.y,fid.x0,fid.R,fid.Q,fid.S,fid.P01,'steadystate','yes');

close all
plotTime(fid.t,fid.x,xfilt,'comp',1:6);
plotFreq(fid.t,fid.x,xfilt,'xlim',[0 50],'comp',1:6);

plotTime(fid.t,fid.p,pfilt);
plotFreq(fid.t,fid.p,pfilt,'xlim',[0 50]);

tilefigs;


