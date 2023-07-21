%% TF of JIS (run Verify_alg first)

[f_axis,Gy]=fft_function(sim.y,mod.dt);
Gy=permute(Gy,[1 3 2]);

omega_axis=f_axis*2*pi;

[Hpy Hx0y Hx1y]=JIS_tf(mod.A,mod.B,mod.G,mod.J,mod.Q,mod.R,mod.S,mod.dt,omega_axis);

Gp=mtimes3(Hpy,Gy);
Gx_filt=mtimes3(Hx0y,Gy);
Gx_pred=mtimes3(Hx1y,Gy);

[~,p_check]=ifft_function(permute(Gp,[1 3 2]),1/mod.dt);
[~,x_filt_check]=ifft_function(permute(Gx_filt,[1 3 2]),1/mod.dt);
[~,x_pred_check]=ifft_function(permute(Gx_pred,[1 3 2]),1/mod.dt);

close all
plottime(sim.t,sim.p,p_jis,real(p_check));

plottime(sim.t,sim.x,x_jis,real(x_filt_check));

