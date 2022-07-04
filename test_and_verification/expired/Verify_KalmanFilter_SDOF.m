%%

clc
clear all
close all

% rng('Default');

%% Create system

dt=0.05

A=exp(-0.05*dt);
G=2.2;

Q=1.2
R=0.8
S=0.5


Ctot=[Q S ; S.' R]

Ccorr=plotcorr(Ctot);

nt=1e6;

t=[1:nt]*dt;

%% Generate noise

[wv]=mvnrnd(zeros(size(Ctot,1),1),Ctot,nt).';

w=wv(1,:);
v=wv(2,:);

C2=corrcoef2(wv);
C2-Ccorr

%%

close all


x=zeros(size(A,1),nt);
y=zeros(size(G,1),nt);

x0=zeros(1,1);
x(:,1)=x0;


for k=1:nt
    
    x(:,k+1)=A*x(:,k)+w(:,k);
    y(:,k)=G*x(:,k)+v(:,k);
    
end
x=x(:,1:nt);

plotTime(t,x);
% plotTime(t,y,y_true);

tilefigs

%% What is this?

% 
% G_id=G(:,[[1:3] [1:3]+4]);
% A_id=A([[1:3] [1:3]+4],[[1:3] [1:3]+4]);
% 
% P_0_0=[]
% 
% 
% Q_id=Q([[1:3] [1:3]+4],[[1:3] [1:3]+4]);
% 
% R_id=R;
% S_id=S([[1:3] [1:3]+4],:);
% x0_id=zeros(6,1);
% 
% [x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss2]=KalmanFilter(A_id,G_id,Q_id,R_id,S_id,y,x0_id,P_0_0);
% [x_k_N,P_k_N]=RTSSmoother(A_id,x_k_k,x_k_kmin,P_k_k,P_k_kmin);
% 
% close all
% plotTime(t,x([[1:3] [1:3]+4],:),x_k_kmin,x_k_k,x_k_N);
% 
% plotFreq(t,x([[1:3] [1:3]+4],:),x_k_kmin,x_k_k,x_k_N,'xlim',[0 5]);
% 


% G_id2=


%% Kalman filter and RTS smoother

P_0_0=[]

[x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KalmanFilter(A,G,Q,R,S,y,x0,P_0_0,'forcescaling','yes');
[x_k_kb,x_k_kminb,P_k_kb,P_k_kminb,K_k_ssb]=KalmanFilter(A,G,Q,R,S*0,y,x0,P_0_0,'forcescaling','yes');

[x_k_N,P_k_N]=RTSSmoother(A,x_k_k,x_k_kmin,P_k_k,P_k_kmin);
[x_k_Nb,P_k_Nb]=RTSSmoother(A,x_k_kb,x_k_kminb,P_k_kb,P_k_kminb);

% close all
% plotTime(t,x,x_k_kmin,x_k_kminb);
% plotTime(t,x,x_k_k,x_k_kb);
% % 
% plotTime(t,x,x_k_N,x_k_Nb);


R_inv=eye(size(R))/R;
A_star=A-S*R_inv*G;
[x_k_N_test,P_k_N_test]=RTSSmoother(A_star,x_k_k,x_k_kmin,P_k_k,P_k_kmin);


%%

clc
close all

var_filt=std(x-x_k_k,0,2).^2
var_filt_b=std(x-x_k_kb,0,2).^2

var_smooth=std(x-x_k_N,0,2).^2
var_smooth_b=std(x-x_k_Nb,0,2).^2

var_smooth_test=std(x-x_k_N_test,0,2).^2


figure();
hold on
plot(var_filt,'db','Markersize',8);
plot(var_filt_b,'*r','Markersize',8)
plot(var_smooth,'om','Markersize',8);
plot(var_smooth_b,'xg','Markersize',8)
plot(var_smooth_test,'sk','Markersize',8)

title('Empirical variance')
legend({'KF ' 'KF , S=0' 'RTS' 'RTS , S=0'});
ylog;

[P_k_kmin_ss,~,~,info]=idare(A.',G.',Q,R,S); %,'noscaling'
P_k_k_ss=P_k_kmin_ss-P_k_kmin_ss*G.'/(G*P_k_kmin_ss*G.'+R)*G*P_k_kmin_ss;

figure();
hold on
plot(diag(P_k_kmin),'db','Markersize',8);
plot(diag(P_k_kminb),'*r','Markersize',8)
title('P_{k|k-1} variance')
ylog;
plot(diag(P_k_kmin_ss),'xg','Markersize',8);
legend({'KF ' 'KF , S=0' 'dare'});

figure();
hold on
plot(diag(P_k_k),'db','Markersize',8);
plot(diag(P_k_kb),'*r','Markersize',8)
title('P_{k|k} variance')
ylog;
plot(diag(P_k_k_ss),'xg','Markersize',8);legend({'KF ' 'KF , S=0' 'dare'});

figure();
hold on;
plot(var_filt./diag(P_k_k),'db','Markersize',8);
plot(var_filt_b./diag(P_k_kb),'*r','Markersize',8)
title('Ratio empirical variance/filter variance (KF)')
ylog;
legend({'KF ' 'KF , S=0' });

figure();
hold on;
plot(var_smooth./diag(P_k_N),'db','Markersize',8);
plot(var_smooth_b./diag(P_k_Nb),'*r','Markersize',8)
plot(var_smooth_test./diag(P_k_N_test),'sk','Markersize',8)

title('Ratio empirical variance/filter variance (RTS)')
ylog;
legend({'RTS ' 'RTS , S=0' });




tilefigs([3 3],'l');

return
%%

clc

e=y-G*x_k_kmin;
eb=y-G*x_k_kminb;

close all

plotTime(t,e,eb);

Omega=(G*P_k_kmin*G.'+R); Omega_inv=eye(size(Omega))/Omega;
Omegab=(G*P_k_kminb*G.'+R); Omegab_inv=eye(size(Omegab))/Omegab;

for k=1:nt
    LL(1,k)=e(:,k).'*Omega_inv*e(:,k);
    LL(2,k)=eb(:,k).'*Omegab_inv*eb(:,k);
end


LL_part2_sum=sum(LL,2)

LL_part1_sum(1,1)=nt*log(det(P_k_kmin));
LL_part1_sum(2,1)=nt*log(det(P_k_kminb));


plotTime(t,LL);
%%

clc
close all

P01=eye(nx)*1e-3;

[x_filt p_filt P_ss Pp_ss M_ss L_ss] = JIS_trunc_ss(A,B,G,J,y,x0,R,Q,S,P01);
[x_filtb p_filtb P_ssb Pp_ssb M_ssb L_ssb] = JIS_trunc_ss(A,B,G,J,y,x0,R,Q,S*0,P01);


[x_filt2 p_filt2 P_ss2 Pp_ss2 M_ss L_ss] = JIS_trunc_ss2(A,B,G,J,y,x0,R,Q,S,P01);
[x_filt2b p_filt2b P_ss2b Pp_ss2b M_ssb L_ssb] = JIS_trunc_ss2(A,B,G,J,y,x0,R,Q,S*0,P01);


close all
plotTime(t,x,x_filt,x_filtb,x_filt2,x_filt2b);xlimall(gcf,[0 1000])
plotTime(t,p,p_filt,p_filtb,p_filt2,p_filt2b);xlimall(gcf,[0 1000])


delta_p=sum((p_filt-p_filt2).^2,2)
delta_pb=sum((p_filtb-p_filt2b).^2,2)

%%

% Note 2020-15-05:

% JIS with S appears to be the same for Maes (2016) and Yan Niu (2011). 

% KF with S appears to be wrong for Maes (2016) (missing effect of S on Kalman gain) and correct for sources in books etc.

