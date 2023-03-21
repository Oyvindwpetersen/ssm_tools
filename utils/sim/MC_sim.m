function x=MC_sim(omega_axis0,S0,t,n,varargin)

%% Simulate with correct time axis (resample included)
%
% Inputs:
% omega_axis0: frequency axis
% S0: spectral density
% t: time axis for simulation
% n: number of simulations
%
% Outputs:
% x: matrix with simulated time series as rows
%
%% Inputs

p=inputParser;

% Parameters for stability
addParameter(p,'seed',[],@isnumeric) 

parse(p,varargin{:})

seed=p.Results.seed;

%% Simulate

% Find maximum omega 
dt_target=diff(t(1:2));
T_target=t(end);

[domegasim,omega_max]=MC_freqaxis(dt_target,T_target);
f_max=omega_max/(2*pi);
dfsim=domegasim/(2*pi);

% Simulation axis and spectrum
f_axis=[dfsim:dfsim:f_max];
omega_axis=2*pi*f_axis;

S_sim(1:size(S0,1),1:size(S0,2),:)=interp1z(omega_axis0,S0,omega_axis,'linear',0);

if any(any(any(isnan(S_sim))))
    error('Cant be NaN');
end

% Simulate by MC

if ~isempty(seed)
    rng(seed);
end

X=MCCholeskyFast(omega_axis0,S0,'Nsim',n,'domegasim',domegasim,'appendomega',omega_max);

% Cut time series
t_sim=X{1};
[~,ind_cut]=min(abs(t(end)*1.01-t_sim));
t_sim=t_sim(1:ind_cut);

n1=size(S0,1);

if n==1 & n1>1
    x_sim=X{2}(:,1:ind_cut);
elseif n1==1
    x_sim=zeros(n,length(t_sim));
    for k=1:n
        x_sim(k,:)=X{k+1}(1:ind_cut);
    end
end

dt_sim=diff(t_sim(1:2));
Fs=1./diff(t(1:2));

Fs_sim=1./diff(t_sim(1:2));

if Fs_sim<Fs
    error('This should be higher, check what is up');
end

% Resample
x=resample(x_sim.',t_sim,Fs).';

x=x(:,1:length(t));

%%