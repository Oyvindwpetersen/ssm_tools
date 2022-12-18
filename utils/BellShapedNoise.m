function x=BellShapedNoise(f_low,f_up,t,n)

%% Simulate noise from bell-shaped spectrum
%
% Inputs:
% f_low: lower frequency bound
% f_up: upper frequency bound
% t: time axis for simulation
% n: number of simulations
%
% Outputs:
% x: matrix with simulated time series as rows
%
%% Check input

if nargin==3
    n=1;
end

if f_low>f_up
   error('f_low must be lower than f_up'); 
end

%% Simulate

% Find maximum omega 
dt_target=diff(t(1:2));
T_target=t(end);

[domegasim,omega_max]=MC_freqaxis(dt_target,T_target);
f_max=omega_max/(2*pi);

% Simulation axis and spectrum
f_axis=[1e-3:1e-3:f_max];
w_axis=f_axis*2*pi;

f_mean=(f_low+f_up)/2;
w_mean=f_mean*2*pi;

% S(w)=exp(-0.5*(w-w_mean).^2/c^2);
w_up=2*pi*f_up;
c=sqrt((-0.5*(w_up-w_mean).^2)/log(0.05));

S(1,1,:)=exp(-0.5*(w_axis-w_mean).^2/c^2);

% Simulate by MC
X=MCCholeskyFast(w_axis,S,'Nsim',n,'domegasim',domegasim);

% Cut time series
t_sim=X{1};
[~,ind_cut]=min(abs(t(end)*1.01-t_sim));
t_sim=t_sim(1:ind_cut);

x_sim=zeros(n,length(t_sim));
for k=1:n
x_sim(k,:)=X{k+1}(1:ind_cut);
x_sim(k,:)=x_sim(k,:)./std(x_sim(k,:));
end

dt_sim=diff(t_sim(1:2));
Fs=1./diff(t(1:2));

% Resample
x=resample(x_sim.',t_sim,Fs).';

x=x(:,1:length(t));


% plotTime(t_sim,x_sim,t,x);