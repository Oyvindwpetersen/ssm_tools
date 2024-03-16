%%

clc
clear all
close all

%%

rng(1);

freq=[0.05 0.2]
omega=diag(freq)*2*pi;
gamma=2*omega*0.02

dt=0.5;

phi=[1 1;-1 2];

Sa=[0 0 ; 0 0].';
Sd=[1 0 ; 0 1  ].';

Sp=[1 0.5].';

[A B G J Ac Bc F]=ssmod_modal(phi,omega,gamma,Sa,Sd,Sp,dt,'force','modal');

lambda=[0.1 0.01];
Fc=-diag(lambda);
Hc=eye(2);
Lc=eye(2);

Ac_aug=[Ac Bc*Hc ; zeros(2,4) Fc];


%% Covariance for LFM

Qcw=diag([3 2]);
Qcw_aug=blkdiag(zeros(size(Ac)),Lc*Qcw*Lc.');

Qdw_aug=cov_c2d(Ac_aug,Qcw_aug,dt);
Qdw_aug_approx=Qcw_aug*dt;

figure(); hold on;
plot(diag(Qdw_aug),'ob','DisplayName','Diagonal of Qw (aug) (exact)');
plot(diag(Qdw_aug_approx),'xr','DisplayName','Diagonal of Qw (aug) (approx)');
legend show
ylog;

plotcorr(Qdw_aug);

%% Covariance for additional white noise excitation

Qce=diag([0.5 0.8]);
Qce_aug=blkdiag(Bc*Qce*Bc.',zeros(2));

Qde_aug=cov_c2d(Ac_aug,Qce_aug,dt);
Qde_aug_approx=Qce_aug*dt;

figure(); hold on;
plot(diag(Qde_aug),'ob','DisplayName','Diagonal of Qw (aug) (exact)');
plot(diag(Qde_aug_approx),'xr','DisplayName','Diagonal of Qw (aug) (approx)');
legend show

plotcorr(Qde_aug);

tilefigs

%%
clc

close all
XTickLabel={'x_1' 'x_2' 'x_3' 'x_4' 's_1' 's_2'};

plotcovmatrix((Qdw_aug),XTickLabel,XTickLabel,'','');
% plotscriptmain('h',6,'w',6,'name','Qw','path',cd,'labelsize',6,'ticksize',6,'legendsize',6,'titlesize',6,'box','on','format',{'jpg'});

plotcovmatrix((Qde_aug),XTickLabel,XTickLabel,'','');
% plotscriptmain('h',6,'w',6,'name','Qe','path',cd,'labelsize',6,'ticksize',6,'legendsize',6,'titlesize',6,'box','on','format',{'jpg'});

tilefigs

%%

Ad_aug=expm(Ac_aug*dt);

%% Effect of dt on 


dt_mat=10.^[-3:0.1:0];

for k=1:length(dt_mat)

    Qdw_aug=cov_c2d(Ac_aug,Qcw_aug,dt_mat(k));

    ratio(k,:)=diag(Qdw_aug)

end

close all
figure();
plot(dt_mat,ratio);
tilefigs;
ylog;
xlog;

%% Check matrix method vs numerical integral method

Qdw_aug1=cov_c2d(Ac_aug,Qcw_aug,dt,'matrix');
Qdw_aug2=cov_c2d(Ac_aug,Qcw_aug,dt,'integral');

Qde_aug1=cov_c2d(Ac_aug,Qce_aug,dt,'matrix');
Qde_aug2=cov_c2d(Ac_aug,Qce_aug,dt,'integral');

ratio_w=Qdw_aug1./Qdw_aug2
ratio_e=Qde_aug1./Qde_aug2

% diff_Q=Qdw_aug1-Qdw_aug2;

XTickLabel={'x_1' 'x_2' 'x_3' 'x_4' 's_1' 's_2'};

close all
plotcovmatrix(ratio_w,XTickLabel,XTickLabel,'','');
% plotscriptmain('h',6,'w',6,'name','Qw_ratio','path',cd,'labelsize',6,'ticksize',6,'legendsize',6,'titlesize',6,'box','on','format',{'jpg'});

plotcovmatrix(ratio_e,XTickLabel,XTickLabel,'','');
% plotscriptmain('h',6,'w',6,'name','Qe_ratio','path',cd,'labelsize',6,'ticksize',6,'legendsize',6,'titlesize',6,'box','on','format',{'jpg'});

tilefigs



