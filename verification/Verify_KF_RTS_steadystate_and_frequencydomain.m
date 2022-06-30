%%
clc
clear all

rng(1);
omega=diag([3 6]);
gamma=2*omega*0.02

dt=0.01;

phi=eye(2);

Sa=[1 0 ; 0 1  ].';
Sd=[0 0 ; 0 0].';

% Sa=[0 0 ; 0 0].';
% Sd=[1 0 ; 0 1  ].';

% Sp=[1 0.5].';

[A B G J Ac Bc F]=ssmod_modal(phi,omega,gamma,Sa,Sd,[],dt);

Tmax=100;
t=[0:dt:Tmax];

p=randn(2,length(t));

x0=zeros(4,1);

[x,y]=statespaceForward(A,B,G,J*0,[],x0,p);

Q=cov((B*p).');
R=diag([4e-6 1e-6])*10;

S=zeros(4,2);

noise=diag(R).^0.5.*randn(size(R,1),length(t));

ynoise=y+noise;

close all
plotTime(t,x);
plotTime(t,y,ynoise);
clear y
%% Estimate

P_0_0=eye(size(Q));
F=A
H=G

[x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KalmanFilter(F,H,Q,R,S,ynoise,x0,P_0_0,'steadystate','yes');
[x_k_N,P_k_N]=RTSSmoother(F,x_k_k,x_k_kmin,P_k_k,P_k_kmin,'steadystate','yes');

close all

plotTime(t,x,x_k_kmin,x_k_k,x_k_N);
plotFreq(t,x,x_k_kmin,x_k_k,x_k_N,'xlim',[0 10]);


figure(); hold on;
plot(diag(P_k_kmin),'-ob');
plot(diag(P_k_k),'-xr');
plot(diag(P_k_N),'-dm');
plot(std(x-x_k_kmin,0,2).^2,'--ob');
plot(std(x-x_k_k,0,2).^2,'--xr');
plot(std(x-x_k_N,0,2).^2,'--dm');
ylog;

legend({'P_k_kmin (filter)' 'P_k_k (filter)' 'P_k_N (filter)' 'P_k_kmin (empirical)' 'P_k_k (empirical)' 'P_k_N (empirical)'})

%% Check that SS and no SS is equal

[x_k_k_noss,x_k_kmin_noss,P_k_k_noss,P_k_kmin_noss]=KalmanFilter(F,H,Q,R,S,ynoise,x0,P_0_0,'steadystate','no');
[x_k_N_noss,P_k_N_noss]=RTSSmoother(F,x_k_k_noss,x_k_kmin_noss,P_k_k_noss,P_k_kmin_noss,'steadystate','no');

close all

plotTime(t,x_k_kmin,x_k_kmin_noss);
plotTime(t,x_k_k,x_k_k_noss);
plotTime(t,x_k_N,x_k_N_noss);

% plotFreq(t,x,x_k_kmin,x_k_k,x_k_N,'xlim',[0 10]);

delta_p=std(x_k_kmin-x_k_kmin_noss,0,2)./std(x_k_kmin,0,2)
delta_f=std(x_k_k-x_k_k_noss,0,2)./std(x_k_k,0,2)
delta_s=std(x_k_N-x_k_N_noss,0,2)./std(x_k_N,0,2)


%% TF of filter and smoother

P_k_kmin_ss=P_k_kmin;
P_k_k_ss=P_k_k;
N_k_ss=P_k_k_ss*F.'/P_k_kmin_ss;

[f_axis,Gy]=fft_function(ynoise,dt);
w_axis=f_axis*2*pi;

% w_axis=[0:0.0001:1]*2*pi;
nx=size(F,1);
ny=size(H,1);

H_yx=zeros(nx*3,ny,length(w_axis));

for k=1:length(w_axis)
    z=exp(1i*w_axis(k)*dt);
    Mat=[
        -eye(nx)+K_k_ss*H eye(nx) zeros(nx) ;
        z*eye(nx) -F zeros(nx) ;
        N_k_ss*z -eye(nx) eye(nx)-N_k_ss*z 
        ];
    
    H_yx(:,:,k)=Mat\[K_k_ss ; zeros(nx,ny) ; zeros(nx,ny) ];

end

close all
plotTransferFunction(w_axis,H_yx,'xlim',[-2 2]*2*pi);
ylimall(gcf,[-2 2])


Hp=H_yx(1:nx,:,:);
Hf=H_yx((nx+1):(nx*2),:,:);
Hs=H_yx((nx*2+1):end,:,:);

% plotTransferFunction(w_axis,Hp);
% plotTransferFunction(w_axis,Hf);
% plotTransferFunction(w_axis,Hs);

plotTransferFunction(w_axis,abs(Hp),w_axis,abs(Hf),w_axis,abs(Hs),'xlim',[-5 5]*2*pi);
plotTransferFunction(w_axis,phasecomplex(Hp),w_axis,phasecomplex(Hf),w_axis,phasecomplex(Hs),'xlim',[-5 5]*2*pi);
linevertical(getSortedAxes(gcf),diag(omega),'--',1,'k');

plotTransferFunction(w_axis,Hp,w_axis,Hf,w_axis,Hs,'xlim',[-5 5]*2*pi);

linevertical(getSortedAxes(gcf),diag(omega),'--',1,'k');

%% Solve in FD

close all

% [f_axis,Gxp]=fft_function(x_k_kmin,dt);
% [f_axis,Gxf]=fft_function(x_k_k,dt);
% [f_axis,Gxs]=fft_function(x_k_N,dt);


Gxp2=zeros(4,length(w_axis));
Gxf2=zeros(4,length(w_axis));
Gxs2=zeros(4,length(w_axis));

for k=1:length(w_axis)
    Gxp2(:,k)=Hp(:,:,k)*Gy(:,k);
    Gxf2(:,k)=Hf(:,:,k)*Gy(:,k);
    Gxs2(:,k)=Hs(:,:,k)*Gy(:,k);
end
    

[tau,x_k_kmin2]=ifft_function(Gxp2,1/dt);
[tau,x_k_k2]=ifft_function(Gxf2,1/dt);
[tau,x_k_N2]=ifft_function(Gxs2,1/dt);


plotTime(t,x_k_kmin,x_k_kmin2);
plotTime(t,x_k_k,x_k_k2);
plotTime(t,x_k_N,x_k_N2);

delta_p=std(x_k_kmin-x_k_kmin2,0,2)
delta_f=std(x_k_k-x_k_k2,0,2)
delta_s=std(x_k_N-x_k_N2,0,2)


