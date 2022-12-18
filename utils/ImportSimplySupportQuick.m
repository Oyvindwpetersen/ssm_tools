function mod=ImportSimplySupportQuick(N_modes)

%%

L=30;
E=210e9;
I=1e-6;
m=50;
N_el=100;

%%

mod=ImportSimplySupport(L,E,I,m,N_el,N_modes)

mod.Xi=eye(size(mod.Omega))*0.02
mod.Gamma=2*mod.Xi.*mod.Omega;


mod.Cg=mod.Mg.*mod.Gamma;


