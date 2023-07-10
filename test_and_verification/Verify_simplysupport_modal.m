%% Verification of frequencies of fe beam vs analytic

clc
clear all
close all

E=210e9
L=5
I=1e-6;
m=50
N_el=100;

n=1:10

f_ref=1/(2*pi)*n.^2*pi^2*sqrt(E*I/(m*L^4));

mod=importbeam(L,E,I,m,N_el);

figure(); hold on; grid on;
title('Natural frequencies');
plot(mod.f(1:10),'ob','DisplayName','Model');
plot(f_ref(1:10),'xr','DisplayName','Analytic');
xlabel('Hz');
legend show


figure(); hold on;

for k=1:10
   plot(mod.phi(1:2:end,k)); 
end

tilefigs