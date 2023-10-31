%%

clc

sim.x0=zeros(mod.nx,1);
[x,y]=ssmod_forward(mod.A,mod.B,mod.G,mod.J,[],sim.x0,sim.p);



[A B O D B_plus B_minus]=ssmod_Galpha(mod.phi,mod.Omega,mod.Gamma,mod.Sa,mod.Sd,mod.Sp,mod.dt,'force','disc');



[x2,y2]=ssmod_forward(A,B,O,D,[],0,sim.p);


close all

plottime(sim.t,y,y2);

Fs=1./mod.dt;


omega_axis=linspace(0.1,0.5*Fs*2*pi,1e5);


H1=ssmod_tf(mod.Ac,mod.Bc,mod.Gc,mod.Jc,omega_axis);

H2=ssmod_tf(mod.A,mod.B,mod.G,mod.J,omega_axis,mod.dt);

H3=ssmod_tf(A,B,O,D,omega_axis,mod.dt);


% plottf(omega_axis,abs(H1),abs(H2),abs(H3),'xlim',[0 100]);

% plottf(omega_axis,angle(H1),angle(H2),angle(H3),'xlim',[0 100]);


plottf(omega_axis,real(H1),real(H2),real(H3),'xlim',[0 100]);

plottf(omega_axis,imag(H1),imag(H2),imag(H3),'xlim',[0 100]);

tilefigs


