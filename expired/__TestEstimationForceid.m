
%%

clc

close all

mod2=struct();

[mod2.A mod2.B mod2.O mod2.D mod2.B_plus mod2.B_minus]=ssmod_Galpha(mod.phi,mod.omega,mod.gamma,mod.Sa,mod.Sd,mod.Sp,dt,'force','modal');

fid2=struct();
fid2.x0=zeros(mod.nm*3,1);

[fid2.xbar fid2.y]=ssmod_forward(mod2.A,mod2.B,mod2.O,mod2.D,[],fid2.x0,fid.p);
fid2.y_clean=fid2.y;

fid2.y=fid2.y_clean+fid.noise;

fid2.x=fid2.xbar+mod2.B_plus*fid.p;

plotTime(fid.t,fid.x,fid2.x(1:mod.nm*2,:));
plotTime(fid.t,fid.y,fid2.y);

plotFreq(fid.t,fid.x,fid2.x(1:mod.nm*2,:),'xlim',[0 10]);
plotFreq(fid.t,fid.y,fid2.y,'xlim',[0 10]);

tilefigs


%%

close all

[xfilt pfilt] = JIS_trunc_ss(mod.A,mod.B,mod.G,mod.J,fid.y,fid.x0,fid.R,fid.Q,fid.S,fid.P01,'trunc','no','plot','yes');

fid2.R=fid.R;
fid2.Q=1e-12*eye(mod.nm*3);
fid2.Q=blockDiagonal(1e-12*eye(mod.nm*2),0*eye(mod.nm*1));
fid2.S=zeros(mod.nm*3,mod.ny);
fid2.P01=fid2.Q;

[xfilt2 pfilt2] = JIS_trunc_ss(mod2.A,mod2.B,mod2.O,mod2.D,fid2.y,fid2.x0,fid2.R,fid2.Q,fid2.S,fid2.P01,'trunc','no','plot','yes','maxsteps',1e6);

[xsmooth2 psmooth2 P_x_ss P_p_ss ]=JIS_smooth(mod2.A,mod2.B,mod2.O,mod2.D,fid2.y,fid2.R,fid2.Q,fid2.S,fid2.x0,fid2.P01,2,'convtol',1e-5,'maxsteps',100e3);


plotTime(fid.t,fid.p,pfilt,pfilt2,'legend',{'True' 'Filt' 'Filt G'});
plotFreq(fid.t,fid.p,pfilt,pfilt2,'xlim',[0 50],'legend',{'True' 'Filt' 'Filt G'});

plotTime(fid.t,fid.p,psmooth,psmooth2,'legend',{'True' 'Filt' 'Filt G'});
plotFreq(fid.t,fid.p,psmooth,psmooth2,'xlim',[0 50],'legend',{'True' 'Filt' 'Filt G'});

% plotTime(fid.t,fid.x(1:6,:),xfilt(1:6,:),xsmooth(1:6,:));
% plotFreq(fid.t,fid.x(1:6,:),xfilt(1:6,:),xsmooth(1:6,:),'xlim',[0 50]);

p_err_SS_JIS=norm(fid.p-pfilt)
p_err_GA_JIS=norm(fid.p-pfilt2)

p_err_SS_SMOOTH=norm(fid.p-psmooth)
p_err_GA_SMOOTH=norm(fid.p-psmooth2)

tilefigs

