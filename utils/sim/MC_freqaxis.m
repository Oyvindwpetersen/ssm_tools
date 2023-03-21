function [domegasim,omega_max]=MC_freqaxis(dt_target,T_target)

%% Find required frequency axis in Monte Carlo simulation with target dt and T
%
% Inputs:
% dt_target: desired maximum time increment
% T_target: desired minimum simulation duration
%
% Outputs:
% domegasim: frequency increment
% omega_max: frequency max
%
%%

domegasim=1e-3;
omega_max=1;
T=2*pi/domegasim;

if T<T_target
    domegasim=2*pi/T_target;
end

converged=false;
while converged==false
    
    % 	omegaaxissim=domegasim:domegasim:omega_max;
    
    N_omega=ceil(omega_max/domegasim);
    NFFT=2^nextpow2(2*N_omega);
    
    % 	t=linspace(0,2*pi/domegasim,NFFT);
    % 	dt=diff(t(1:2));
    
    dt=(2*pi/domegasim)/(NFFT-1);
    
    if dt>dt_target
        omega_max=omega_max+0.1;
    else
        converged=true;
    end
    
end


