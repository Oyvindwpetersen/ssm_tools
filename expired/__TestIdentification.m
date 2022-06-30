%%

clc
clear all
close all

hardangerFEModel
hardangerAeroModel

mod.a_ref={ };
mod.a_o=setdiff(mod.a_o,mod.a_ref);
mod.d_o=mod.a_o;
mod.d_o={};

[mod.Sd,mod.Sa,mod.Sp,]= generalSelection(mod.d_o,mod.a_o,{},basemod.doflabel);

mod.nm=length(mod.includeModes);
mod.np=mod.nm;
mod.ny=size(mod.Sa,1);

%%

mod.M0_g=eye(mod.nm);
mod.C0_g=mod.gamma;
mod.K0_g=mod.omega^2;

mod.M_g=mod.M0_g;
mod.C_g=mod.C0_g;
mod.K_g=mod.K0_g;

[mod.A mod.B mod.G mod.J mod.Ac mod.Bc]=statespaceModelFull(mod.K_g,mod.C_g,mod.M_g,mod.Sa*mod.phi,mod.Sd*mod.phi,eye(mod.nm),mod.dt);

[~,interval.phi_eff{k},interval.omega_eff(:,k),interval.xi_eff(:,k)]=eigA(mod.A,mod.dt);

lambda=1;
sigma_p=100;
sigma_w=sqrt(2*lambda*sigma_p.^2);
sigma_e=1;

lambda_mat_glob=lambda*ones(1,mod.nm)';
sigma_mat_glob=sigma_w*ones(1,mod.nm)';

[mod.Fc,mod.Lc,mod.Hc]=buildMaternModelGlobal(lambda_mat_glob,sigma_mat_glob,0,[],[],[])

[mod.Fac,mod.Bac,mod.Hac,mod.Jac,mod.Fad,mod.Bad,mod.Had,mod.Jad]=statespaceModelLatentForce(mod.Ac,mod.Bc,mod.G,mod.J,mod.Fc,mod.Hc,mod.Lc,NaN,NaN,mod.dt)

mod.Qwad=covarianceContToDisc(mod.Fac,blockDiagonal(zeros(mod.nm*2),mod.Lc*sigma_w^2*eye(mod.nm)*mod.Lc.'),mod.dt);

mod.Qead=covarianceContToDisc(mod.Fac,blockDiagonal(mod.Bc*sigma_e^2*eye(mod.nm)*mod.Bc.',zeros(mod.nm)),mod.dt);

%%

clc
close all

fid.t=[0:mod.dt:1800];

wk=mvnrnd(zeros(size(mod.Qwad,1),1),mod.Qwad,length(fid.t)).';
epsk=mvnrnd(zeros(size(mod.Qead,1),1),mod.Qead,length(fid.t)).';

mod.R=eye(mod.ny)*1e-6;

vk=mvnrnd(zeros(size(mod.R,1),1),mod.R,length(fid.t)).';

% xa_true=zeros(mod.nm*3,length(fid.t));
% for k=1:length(fid.t)-1
%     xa_true(:,k+1)=mod.Fad*xa_true(:,k)+wk(:,k)+epsk(:,k)*0;
% end
% p_true=xa_true(mod.nm*2+1:end,:);


p_true=zeros(mod.np,length(fid.t));
for k=1:mod.np
    p_true(k,:)=sawtooth(2*pi*50*fid.t);
end

for k=1:length(fid.t)-1
    xa_true(:,k+1)=mod.Fad*xa_true(:,k)+wk(:,k)+epsk(:,k)*0;
end
p_true=xa_true(mod.nm*2+1:end,:);


% plotTime(fid.t,xa_true(1:mod.nm*2,:),'comp',[1:10]);

% plotTime(fid.t,p_true);

fid.y_true=mod.Had*xa_true;
fid.y=fid.y_true+vk;

close all

% plotTime(fid.t,mod.phi_acc*xa_true(1:mod.nm,:));
% plotFreq(fid.t,mod.phi_acc*xa_true(1:mod.nm,:));
% 
% plotTime(fid.t,fid.y);
% plotFreq(fid.t,fid.y);

tilefigs


%%

close all

mod.Qad=mod.Qwad;

Q_test=blockDiagonal(eye(mod.nm*2),zeros(mod.nm))*100000;

fid.x0=zeros(mod.nm*3,1);
fid.P_0_0=zeros(mod.nm*3);

mod.S=zeros(mod.nm*3,mod.ny);

[xa_hat]=KalmanFilter(mod.Fad,mod.Had,mod.Qwad+Q_test,mod.R,mod.S,fid.y,fid.x0,fid.P_0_0);

p_hat=xa_hat(mod.nm*2+1:end,:);

plotTime(fid.t,xa_true,xa_hat,'comp',[1:mod.nm]);

plotTime(fid.t,p_true,p_hat);
plotFreq(fid.t,p_true,p_hat);

tilefigs

