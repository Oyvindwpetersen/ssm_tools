%%

clc
clear all
close all

L=100
E=210e9
I=0.05
m=8000
N_el=100

[Kg,Mg,phi,phi_full,K,M,K_red,M_red,S_red]=SimplySupportModal(L,E,I,m,N_el);

figure();
plot(phi_full(1:2:end,1:5));

wn=sqrt(full(diag(Kg)./diag(Mg)));

fn=wn/(2*pi)

%%

doflabel=getLabel({'U' 'UR'},[1:(N_el+1)]);

u_label=getLabel({'U'},[1:(N_el+1)]);
p_label=getLabel({'U'},[20]);

[Sd,~,Sp]=DofSelection(u_label,{},p_label,doflabel);

% Sp=Sp*ones(N_el+1,1);
% Sp=Sp*ones(N_el/2,1);

ps=100

rs_all=[]
for k=1:100
    
zs=Kg(1:k,1:k)\phi_full(:,1:k).'*Sp*ps;

rs=Sd*phi_full(:,1:k)*zs

rs_all(:,k)=rs;

end

close all

figure(); hold on;
plot(rs)

figure(); hold on;
plot(rs_all(:,1:10));

%%

doflabel=getLabel({'U' 'UR'},[1:(N_el+1)]);

mod1=struct();
mod1.nm=100
mod1.omega=Kg(1:mod1.nm,1:mod1.nm).^0.5;
mod1.gamma=2*mod1.omega*0.01;
mod1.phi=phi_full(:,1:mod1.nm);

mod2=struct();
mod2.nm=10
mod2.omega=Kg(1:mod2.nm,1:mod2.nm).^0.5;
mod2.gamma=2*mod2.omega*0.01;
mod2.phi=phi_full(:,1:mod2.nm);

close all

dt=0.001;
d_cell={'20_U'} % '30_U' '40_U' '50_U' '60_U'
a_cell={'20_U' } %'30_U' '40_U' '50_U' '60_U'
p_cell={'20_U'}
[Sd,Sa,Sp]=DofSelection(d_cell,a_cell,p_cell,doflabel);

[mod1.A mod1.B mod1.G mod1.J mod1.Ac mod1.Bc]=ssmod_modal(mod1.phi,mod1.omega,mod1.gamma,Sa,Sd,Sp,dt,'force','disc');
[mod2.A mod2.B mod2.G mod2.J mod2.Ac mod2.Bc]=ssmod_modal(mod2.phi,mod2.omega,mod2.gamma,Sa,Sd,Sp,dt,'force','disc');

Tmax=30;
t=[0:dt:Tmax];

p=trianglepulse(2.4,2.5,2.6,t)*1e3;

plotTime(t,p);

mod1.Q=1e-14*eye(size(mod1.A));
mod1.R=1e-12*eye(size(mod1.G,1));
mod1.S=zeros(size(mod1.Q,1),size(mod1.R,1));

[w,v]=cov_noisegen(mod1.Q,mod1.R,mod1.S,t);

mod1.x0=zeros(size(mod1.A,1),1);

[x1,y1_clean]=ssmod_forward(mod1.A,mod1.B,mod1.G,mod1.J,[],mod1.x0,p);

y1=y1_clean+v;

plotTime(t,y1);
plotFreq(t,y1,'xlim',[0 100]);

return
%%


% mod1.P_0_1=eye(size(mod1.Q));
% mod1.Q=diag([ 1e-3*ones(1,20) 1e-3*ones(1,180) ]);

% mod1.Q=blkdiag(mod1.Q,mod1.Q);
mod1.P_0_1=eye(size(mod1.Q));

[fid1.x_jis fid1.p_jis]=JIS_trunc_ss(sparse(mod1.A),sparse(mod1.B),mod1.G,mod1.J,y1,mod1.x0,mod1.R,mod1.Q,sparse(mod1.S),mod1.P_0_1,'trunc','no','dispconv',true,'maxsteps',200e3,'scale','yes');

% [fid1.x_smooth fid1.p_smooth Px_smooth_ss Pp_smooth_ss ]=JIS_smooth(mod1.A,mod1.B,mod1.G,mod1.J,y1,mod1.R,mod1.Q,mod1.S,mod1.x0,mod1.P_0_1,3,'convtol',1e-6,'maxsteps',200e3);

% % 
close all
% % 
% plotTime(t,x1,fid1.x_jis,'comp',[1:10]);
% plotTime(t,p,fid1.p_jis);


% return

mod2.Q=1e-6*eye(size(mod2.A));
mod2.R=mod1.R;
mod2.S=zeros(size(mod2.Q,1),size(mod2.R,1));
mod2.P_0_1=eye(size(mod2.Q));

mod2.x0=zeros(size(mod2.A,1),1);

[fid2.x_jis fid2.p_jis]=JIS_trunc_ss(mod2.A,mod2.B,mod2.G,mod2.J,y1,mod2.x0,mod2.R,mod2.Q,mod2.S,mod2.P_0_1,'maxsteps',200e3,'scale','yes');


% [fid1.x_smooth fid1.p_smooth Px_smooth_ss Pp_smooth_ss ]=JIS_smooth(mod2.A,mod2.B,mod2.G,mod2.J,y1,mod2.R,mod2.Q,mod2.S,mod2.x0,mod2.P_0_1,3,'convtol',1e-8);

plotTime(t,p,fid1.p_jis,fid2.p_jis);
plotTime(t,p,fid2.p_jis);


return
% [Sd,~,~]=DofSelection(u_label,{},{},doflabel);
% [~,~,Sps]=DofSelection({},{},ps_label,doflabel);
% [~,~,Spd]=DofSelection({},{},pd_label,doflabel);






% u_label=getLabel({'U'},[1:(N_el+1)]);
% ps_label=getLabel({'U'},[20]);
% pd_label=getLabel({'U'},[1:(N_el+1)]);
% 
% [Sd,~,~]=DofSelection(u_label,{},{},doflabel);
% [~,~,Sps]=DofSelection({},{},ps_label,doflabel);
% [~,~,Spd]=DofSelection({},{},pd_label,doflabel);


% ps0=100
% 
% sim1.zs=mod1.omega.^2\mod1.phi.'*Sps*ps0;
% sim1.us=Sd*mod1.phi*sim1.zs;
% 
% sim2.zs=mod2.omega.^2\mod2.phi.'*Sps*ps0;
% sim2.us=Sd*mod2.phi*sim2.zs;
% 
% figure(); hold on;
% plot(sim1.us,'-b')
% plot(sim2.us,'--r')


%%

dt=0.05;
Tmax=300;

t=[0:dt:Tmax];

ps=BandLimitedWhiteNoise(1e-3,1e-2,t,1)+10;
pd=BandLimitedWhiteNoise(1e-2,1e0,t,N_el+1);

T=exp(-abs([0:(N_el)]-[0:(N_el)].')/10); T=forcesym(T);
[V,D]=eig(T); T2=V*D.^0.5;
pd2=T2*pd;

plotcorr(corrcoef(pd2.'))

plotTime(t,ps,pd(1,:),ps+pd(1,:))

%%


