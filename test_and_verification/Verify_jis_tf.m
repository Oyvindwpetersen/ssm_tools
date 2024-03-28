%% TF of JIS (run Verify_algorithms first)


[x_jis p_jis Px_jis_ss Pp_jis_ss] = JIS_ss(mod.A,mod.B,mod.G,mod.J,sim.y,sim.x0,mod.Q,mod.R,mod.S,mod.P01,'trunc',false);



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


% Cut to remove initial effects

plotfreq(sim.t(:,1000:end),sim.x(:,1000:end),x_jis(:,1000:end),real(x_filt_check(:,1000:end)));
plotfreq(sim.t(:,1000:end),sim.p(:,1000:end),p_jis(:,1000:end),real(p_check(:,1000:end)));


