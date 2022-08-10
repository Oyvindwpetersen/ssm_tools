%%



[fid1.x_jis fid1.p_jis]=JIS_trunc_ss(sparse(mod1.A),sparse(mod1.B),mod1.G,mod1.J,y1,mod1.x0,mod1.R,mod1.Q,sparse(mod1.S),mod1.P_0_1,'trunc','yes','maxsteps',10000);
% [fid1.x_smooth fid1.p_smooth Px_smooth_ss Pp_smooth_ss ]=JIS_smooth(mod1.A,mod1.B,mod1.G,mod1.J,y1,mod1.R,mod1.Q,mod1.S,mod1.x0,mod1.P_0_1,3,'convtol',1e-6,'maxsteps',200e3);


As=sparse(mod1.A);
Bs=sparse(mod1.B);

tic 
for k=1:1000
    c=[As' ; Bs'];
end
toc