%%

clc
clear all
close all

omega_axis=linspace(0,5,1e3);

m=1;
d=0.6;
a=-0.15;

c=0.02;

k=4;

H_str(1,1,:)=1./(-omega_axis.^2*m+i*omega_axis*c+k);

% plottf(omega_axis,abs(H_str));

[~,~,~,~, Ac0]=ssmod_full(k,c,m,1,0,1,0.01);

[lambda0,phi0,omega0,xi0]=eigA(Ac0);

G(1,1,:)=a./(i*omega_axis+d)

H_tot(1,1,:)=1./(-omega_axis.^2*m+i*omega_axis.*(c+a./(i*omega_axis+d))+k);


% plottf(omega_axis,abs(H_str),abs(H_tot),abs(G));

p1=[m , m*d+c , a+c*d+k , k*d]

lambda1=roots(p1)
omega1=abs(lambda1)
xi1=-real(lambda1)./abs(lambda1)


lam_guess=0;
for ind=1:10

p2=[m-a/(abs(lam_guess).^2+d^2) , c+a*d/(abs(lam_guess).^2+d^2) , k];

r=roots(p2);
lam_guess=r(1);

lam_guess_save(ind,1)=lam_guess;

end

lambda2=lam_guess;
omega2=abs(lambda2)
xi2=-real(lambda2)./abs(lambda2)



ratio_om=omega2./omega1(1)
ratio_xi=xi2./xi1(1)
