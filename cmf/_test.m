%%

clc
clear all
close all

w_axis=logspace(-6,3,2000);

sigma=2.2
Au=15*0.1
Ku=10*0.1
U=15
z=60

Su=KaimalSpectrum(w_axis,Au,U,sigma,z);

dt=0.1/10

j_lag=2.^[0:4];
l_lag=2.^[0:4];

S_target(1,1,:)=Su;

n=11
L=2.8
x_axis=linspace(0,L,n);
clear S_target
for k1=1:n
for k2=1:n
    dx=abs(diff(x_axis([k1 k2])));
    S_target(k1,k2,:)=Su.*exp(-w_axis*dx*Ku);
end
end

[A_mat,B,R_target,tau_target]=AR_fit(w_axis,S_target,j_lag,l_lag,dt);

lag_calc=[0:2000]; tau_target=lag_calc*dt;
R_target=zeros(size(S_target,1),size(S_target,1),length(lag_calc));

for k=length(lag_calc):-1:1
        cos_vec(1,1,:)=cos(w_axis*dt*lag_calc(k));
        cos_mat=repmat(cos_vec,size(S_target,1),size(S_target,1),1);
        R_target(:,:,k)=simps(w_axis,S_target.*cos_mat,3);
end


[S_AR,H]=AR_PSD(w_axis,A_mat,B,j_lag,dt);

plotSpectrum(w_axis,S_target,w_axis,S_AR*2,'comp',[1,5,11],'xlim',[0 2*2*pi]);


l_lag2=[0:2000];

[R_AR]=AR_CF(A_mat,B,j_lag,l_lag2,5000);

tau_AR=l_lag2*dt;

plotAutocorr(tau_target,R_target,tau_AR,R_AR,'comp',[1],'xlim',[0 100]);


% plotSpectrum(w_axis,S2coh(S_target),w_axis,S2coh(S_AR),'comp',[1 2 5 10 ],'xlim',[0 5*2*pi]);

% plotSpectrum(w_axis,S,w_axis,abs(S_AR),'comp',[1 5 11],'xlim',[0 5*2*pi]);


tilefigs([],'l');
tilefigs
