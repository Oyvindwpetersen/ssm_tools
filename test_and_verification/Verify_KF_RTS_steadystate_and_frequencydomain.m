%% Verification of KF/RTS

clc
clear all
close all

% rng('Default');

%% Create system

mod=importbeamquick(2);

mod.dt=0.05;

%% Disc force

mod.a_cell={'10_U' '20_U' '30_U' '40_U' '50_U' '60_U' '70_U' '80_U' '90_U'}
mod.d_cell=mod.a_cell;
mod.p_cell={'30_U'}
[mod.Sd,mod.Sa,mod.Sp]=DofSelection(mod.d_cell,mod.a_cell,mod.p_cell,mod.doflabel);

[mod.A mod.B mod.G mod.J mod.Ac mod.Bc]=ssmod_modal(mod.phi,mod.Omega,mod.Gamma,mod.Sa,mod.Sd,mod.Sp,mod.dt,'force','disc');

mod.ny=size(mod.Sa,1); mod.nx=size(mod.A,1); mod.np=size(mod.Sp,2);

%%

mod.Q=eye(mod.nx)*[1e-2]^2;
mod.R=eye(mod.ny)*[1e-4]^2;
mod.S=zeros(mod.nx,mod.ny);
% mod.S=[ones(mod.nx/2,mod.ny)*0.1 ; ones(mod.nx/2,mod.ny)*-0.1].*diag(mod.Q).^0.5.*(diag(mod.R).').^0.5;

Ctot=[mod.Q mod.S ; mod.S.' mod.R]
Ctot=plotcorr([mod.Q mod.S ; mod.S.' mod.R]);

sim.nt=1e4+1;
sim.t=[1:sim.nt]*mod.dt;


%% Generate noise

[sim.w,sim.v]=cov_noisegen(mod.Q,mod.R,mod.S,sim.t);

close all

sim.x0=zeros(mod.nx,1);
[sim.x,sim.y]=ssmod_forward_stoch(mod.A,[],mod.G,[],[],sim.x0,[],sim.w,sim.v);

% plottime(sim.t,sim.x);
% plottime(sim.t,sim.y);

tilefigs

%% Estimate

P_0_0=eye(size(mod.Q));

[x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KalmanFilter(mod.A,mod.G,mod.Q,mod.R,mod.S,sim.y,sim.x0,P_0_0,'steadystate','yes');
[x_k_N,P_k_N]=RTSSmoother(mod.A,x_k_k,x_k_kmin,P_k_k,P_k_kmin,'steadystate','yes');

close all

plottime(sim.t,sim.x,x_k_kmin,x_k_k,x_k_N);
plotfreq(sim.t,sim.x,x_k_kmin,x_k_k,x_k_N,'xlim',[0 10]);


figure(); hold on;
plot(diag(P_k_kmin),'-ob');
plot(diag(P_k_k),'-or');
plot(diag(P_k_N),'-om');
plot(std(sim.x-x_k_kmin,0,2).^2,'--xb');
plot(std(sim.x-x_k_k,0,2).^2,'--xr');
plot(std(sim.x-x_k_N,0,2).^2,'--xm');
ylog;

legend({'P_k_kmin (filter)' 'P_k_k (filter)' 'P_k_N (filter)' 'P_k_kmin (empirical)' 'P_k_k (empirical)' 'P_k_N (empirical)'})

%% Check that SS and no SS is equal

[x_k_k_noss,x_k_kmin_noss,P_k_k_noss,P_k_kmin_noss]=KalmanFilter(mod.A,mod.G,mod.Q,mod.R,mod.S,sim.y,sim.x0,P_0_0,'steadystate','np');
[x_k_N_noss,P_k_N_noss]=RTSSmoother(mod.A,x_k_k_noss,x_k_kmin_noss,P_k_k_noss,P_k_kmin_noss,'steadystate','no');

close all

plottime(sim.t,x_k_kmin,x_k_kmin_noss);
plottime(sim.t,x_k_k,x_k_k_noss);
plottime(sim.t,x_k_N,x_k_N_noss);

% plotfreq(t,x,x_k_kmin,x_k_k,x_k_N,'xlim',[0 10]);

delta_p=std(x_k_kmin-x_k_kmin_noss,0,2)./std(x_k_kmin,0,2)
delta_f=std(x_k_k-x_k_k_noss,0,2)./std(x_k_k,0,2)
delta_s=std(x_k_N-x_k_N_noss,0,2)./std(x_k_N,0,2)


%% TF of filter and smoother

P_k_kmin_ss=P_k_kmin;
P_k_k_ss=P_k_k;
N_k_ss=P_k_k_ss*mod.A.'/P_k_kmin_ss;

[f_axis,Gy]=fft_function(sim.y,mod.dt);
w_axis=f_axis*2*pi;

% w_axis=[0:0.0001:1]*2*pi;
nx=size(mod.A,1);
ny=size(mod.G,1);

H_yx=zeros(nx*3,ny,length(w_axis));

for k=1:length(w_axis)
    z=exp(1i*w_axis(k)*mod.dt);
    Mat=[
        -eye(nx)+K_k_ss*mod.G eye(nx) zeros(nx) ;
        z*eye(nx) -mod.A zeros(nx) ;
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

Gxp2=zeros(mod.nx,length(w_axis));
Gxf2=zeros(mod.nx,length(w_axis));
Gxs2=zeros(mod.nx,length(w_axis));

for k=1:length(w_axis)
    Gxp2(:,k)=Hp(:,:,k)*Gy(:,k);
    Gxf2(:,k)=Hf(:,:,k)*Gy(:,k);
    Gxs2(:,k)=Hs(:,:,k)*Gy(:,k);
end

[tau,x_k_kmin2]=ifft_function(Gxp2,1/mod.dt);
[tau,x_k_k2]=ifft_function(Gxf2,1/mod.dt);
[tau,x_k_N2]=ifft_function(Gxs2,1/mod.dt);

plottime(sim.t,x_k_kmin,x_k_kmin2);
plottime(sim.t,x_k_k,x_k_k2);
plottime(sim.t,x_k_N,x_k_N2);

delta_p=std(x_k_kmin-x_k_kmin2,0,2)
delta_f=std(x_k_k-x_k_k2,0,2)
delta_s=std(x_k_N-x_k_N2,0,2)


