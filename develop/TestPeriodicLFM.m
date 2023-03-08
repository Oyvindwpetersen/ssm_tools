%%

clc
clear all
close all

lambda=0.01

w0=0.1*2*pi;

dt=0.05;
T=6000

t=[0:dt:T];

sigma_p=3.4
ns=3;

[Fc,Lc,Hc,Qc,sigma_w]=ssmod_periodicdecay(w0,lambda,sigma_p,'matern',ns);


[Fd,~,Hd,~]=ssmod_c2d(Fc,Lc,Hc,0,dt);

Bd=eye(size(Fd,1));

Qd=cov_c2d(Fc,Lc*Qc*Lc.',dt);



%%

close all

x0=zeros(size(Qd,1),1);

w=mvnrnd(zeros(size(Qd,1),1),Qd,length(t)).';

[x_sim,y_sim]=ssmod_forward(Fd,Bd,Hd,zeros(1,size(Qd,1)),[],x0,w);

% plotTime(t,x_sim);

plotTime(t,y_sim);
plotFreq(t,y_sim);

tilefigs


[R,tau]=xcorrSignal(y_sim,dt,100000);
% 
% 
% Theory
R_theory(1,1,:)=sigma_w^2*exp(-lambda*tau)/(2*lambda).*cos(w0*tau);
% 
plotAutocorr(tau,R,tau,R_theory,'xlim',[0 200]);
%
% 

%%

clc
clear all
close all


tau=0:0.01:100;

for j=0:100;
    
    k(j+1,:)=1/factorial(j)*cos(tau).^j;
    
end


figure(); hold on; grid on;
plot(tau,sum(k,1));
plot(tau,sum(k(2:end,:),1));


%%




% Fq=-lambda;
% Lq=1;
% Hq=1;
% [Fq,Lq,Hq]=ssmod_squaredexp(1/lambda,1,3,[]);
% 
% Fp=[0 -w0*1 ; w0*1 0];
% Lp=eye(2);
% Hp=[1 0];

% [Fp_d,Lp_d

% q=size(Fq,1);
% 
% Fc=kron(Fq,eye(2))+kron(eye(q),Fp);
% 
% Lc=kron(Lq,Lp);
% 

% qj=1
% Sigma_w=3.3
% 
% Qc=kron(Sigma_w^2,qj^2*eye(2));
% 
% Hc=kron(Hq,Hp);

% 
% 
% ssmod_perioddecay(w0,lambda,sigma_p,kernel,ns)
% 
