%%

clc
clear all
close all

lambda=0.1
Ac=-lambda;
Bc=1;
Gc=1;
Jc=1;
[A,B,G,J]=ssmod_c2d(Ac,Bc,Gc,Jc,dt);
R=0.001;
Q=1e-12;
S=0;








dt=0.01;
nt=1e6;
t=[0:nt]*dt;

p=BandLimitedWhiteNoise(1,1.5,t,1);




w=randn(size(t))*sqrt(Q);
v=randn(size(t))*sqrt(R);

x0=0;
P_0_1=Q;

[x,y]=ssmod_forward_stoch(A,B,G,J,[],x0,p,w,v);

plotTime(t,y);

clc
close all

[x_smooth p_smooth Px_smooth_ss Pp_smooth_ss] = JIS_smooth(A,B,G,J,y,R,Q,S,x0,P_0_1,10,'convtol',1e-8);

[x_smooth2 p_smooth2 Px_smooth_ss2 Pp_smooth_ss2] = JIS_smooth_testnew(A,B,G,J,y,R,Q,S,x0,P_0_1,10,'convtol',1e-8);

var_state_smooth=std(x-x_smooth,0,2).^2
var_force_smooth=std(p-p_smooth,0,2).^2

% var_state_smooth2=std(x-x_smooth2,0,2).^2
% var_force_smooth2=std(p-p_smooth2,0,2).^2

var_state_smooth=std(x(:,1:end-L)-x_smooth(:,1:end-L),0,2).^2
var_force_smooth=std(p(:,1:end-L)-p_smooth(:,1:end-L),0,2).^2

figure(); hold on; grid on; ylog;
title('Ratio variance empirical/predicted'); legend show; 
plot(var_state_smooth./diag(Px_smooth_ss),'db','Markersize',8,'DisplayName','Smoother (state)');
plot(var_force_smooth./diag(Pp_smooth_ss),'*r','Markersize',8,'DisplayName','Smoother (force)');


figure(); hold on; grid on; ylog;
title('Ratio variance empirical/predicted'); legend show; 
plot(var_state_smooth2./diag(Px_smooth_ss2),'db','Markersize',8,'DisplayName','Smoother (state)');
plot(var_force_smooth2./diag(Pp_smooth_ss2),'*r','Markersize',8,'DisplayName','Smoother (force)');

tilefigs([3 4],'l');

