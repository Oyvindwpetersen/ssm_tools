%%


clc


mod.Hc=1;
mod.Lc=1;

lambda=0.01
sigma_w=2

mod.Fc=-lambda

mod.Fac=[mod.Ac mod.Bc*mod.Hc ; zeros(1,size(mod.Ac,2)) mod.Fc];

mod.Hac=[mod.G mod.J*mod.Hc];

Qwc=sigma_w

fid.Qad=blkdiag(zeros(size(mod.Ac)),Qwc*dt);

[mod.Fad,~,mod.Had]=statespaceContToDisc(mod.Fac,[],mod.Hac,[],dt)

varR=var(fid.y,0,2)*0.01^2;
fid.R=diag(varR);
fid.Q=eye(mod.nm*2)*1e-14;


[fid.x fid.y]=statespaceForward(mod.A,mod.B,mod.G,mod.J,mod.F*0,fid.x0,fid.p);

fid.y_clean=fid.y;
rng(0); fid.y=fid.y_clean+randn(mod.nd,length(fid.t)).*varR.^0.5;

addPathTemp2('C:\Cloud\OD_OWP\Work\Hardanger\GPLFM\latentforcemodel');

fid.Sa=zeros(size(fid.Qad,1),size(fid.R,1));

fid.xa0=zeros(size(mod.Fad,1),1);
fid.Pa_0_0=zeros(size(mod.Fad,1));

[xa_est,~,~,~,~]=KalmanFilter(mod.Fad,mod.Had,fid.Qad,fid.R,fid.Sa,fid.y,fid.xa0,fid.Pa_0_0);

x_filt2=xa_est(1:end-1,:);
p_filt2=xa_est(end,:);


fid.P01=fid.Q;
fid.x0=zeros(mod.nm*2,1);
fid.S=zeros(size(fid.Q,1),size(fid.R,1));

[x_filt p_filt] = JIS_trunc(mod.A,mod.B,mod.G,mod.J,fid.y,fid.x0,fid.R,fid.Q,fid.S,fid.P01,'steadystate','yes');


close all


plotopt=struct;
plotopt.LineStyleSet={'-' '--' '--' ':'}


% 
plotTime(fid.t,fid.x,x_filt,x_filt2,'comp',[1:6],plotopt);
plotTime(fid.t,fid.p,p_filt,p_filt2,plotopt);

tilefigs([3 3 ])

norm(p_filt-fid.p)
norm(p_filt2-fid.p)


%%


lambda_mat=10.^[-4:3];
sigma_p_mat=10.^[-5:2];


clear p_norm
for ind1=1:length(lambda_mat)
for ind2=1:length(sigma_p_mat)
    
    
    lambda=lambda_mat(ind1);
    sigma_w=sigma_p_mat(ind2)*sqrt(2*lambda);
    
    mod.Fc=-lambda

    mod.Fac=[mod.Ac mod.Bc*mod.Hc ; zeros(1,size(mod.Ac,2)) mod.Fc];

    mod.Hac=[mod.G mod.J*mod.Hc];

    Qwc=sigma_w

    fid.Qad=blkdiag(zeros(size(mod.Ac)),Qwc*dt);

    [mod.Fad,~,mod.Had]=statespaceContToDisc(mod.Fac,[],mod.Hac,[],dt)

    varR=var(fid.y,0,2)*0.1^2;
    fid.R=diag(varR);
    fid.Q=eye(mod.nm*2)*1e-14;

    [fid.x fid.y]=statespaceForward(mod.A,mod.B,mod.G,mod.J,mod.F*0,fid.x0,fid.p);

    fid.y_clean=fid.y;
    rng(0); fid.y=fid.y_clean+randn(mod.nd,length(fid.t)).*varR.^0.5;

    addPathTemp2('C:\Cloud\OD_OWP\Work\Hardanger\GPLFM\latentforcemodel');

    fid.S=zeros(size(fid.Qad,1),size(fid.R,1));

    fid.xa0=zeros(size(mod.Fad,1),1);
    fid.Pa_0_0=zeros(size(mod.Fad,1));

    [xa_est,~,~,~,~]=KalmanFilter(mod.Fad,mod.Had,fid.Qad,fid.R,fid.S,fid.y,fid.xa0,fid.Pa_0_0);

    x_filt2=xa_est(1:end-1,:);
    p_filt2=xa_est(end,:);
    
    
    
    
    
    
    p_norm(ind1,ind2)=norm(p_filt2-fid.p);
    
    
    
    
    
    
end
end

close all

figure();
[l_grid,p_grid]=meshgrid(lambda_mat,sigma_p_mat);
surf(l_grid,p_grid,p_norm.'); ylog; xlog; zlog;
xlabel('lam');
ylabel('p');
