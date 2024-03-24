function [Fc_out,Lc_out,Hc_out,lambda_hat,sigma_w_hat,sigma_p_hat,S_par_out,R_par_out,ind]=fitmatern_mom(tau_axis,R_target,omega_axis,S_target,varargin)

%%

p=inputParser;
addParameter(p,'plot','no',@ischar)
addParameter(p,'h',10,@isnumeric)
addParameter(p,'maxtau',100,@isnumeric)
addParameter(p,'omegan',[],@isnumeric)
addParameter(p,'scaletopsd',false,@islogical)

parse(p,varargin{:});

doPlot=p.Results.plot;
h=p.Results.h;
maxtau=p.Results.maxtau;
omegan=p.Results.omegan;
scaletopsd=p.Results.scaletopsd;

%% Overdetermined system

if ndims(R_target)==2
    R_target=diag3d(R_target);
end

p_mat=[0];

tau_m=[0:h:maxtau];
ind=indexargmin(tau_axis,tau_m);

diff_tau=tau_axis(ind)-tau_m;
if any(diff_tau~=0)
    error('Tau axis of chosen points must be on the input tau axis');
end

for k=1:size(R_target,1)
    
    R_temp=squeeze(R_target(k,k,ind));

    clear Gamma_m gamma_m
    Gamma_m(:,1)=R_temp(1:end-1);    
    gamma_m(:,1)=R_temp(2:end);

    phi_hat(k,1)=Gamma_m\gamma_m;
    
    if phi_hat(k,1)<0
         phi_hat(k,1)=- phi_hat(k,1);
         warning(['Sign of phi in mom reversed for mode=' num2str(k)]);
    end    
    
    lambda_hat(k,1)=-1./h*log(phi_hat(k));
    
    sigma_p_hat(k,1)=sqrt(R_temp(1));
    sigma_w_hat(k,1)=sqrt(2*lambda_hat(k)*sigma_p_hat(k).^2);

    
end

if scaletopsd
    for k=1:size(R_target,1)

        % Divide by 2 since S_target is a one-sided spectrum
        S_target_at_freq=interp1z(omega_axis,S_target(k,k,:),omegan(k)) /2;
        
        % Two sided Matern: 
        sigma_w_hat_adjusted(k,1)=sqrt( S_target_at_freq*2*pi*(lambda_hat(k,1).^2+omegan(k).^2) );
        sigma_p_hat_adjusted(k,1)=sqrt( sigma_w_hat_adjusted(k,1).^2/(2*lambda_hat(k,1)) );
        
        adjustment_factor(k,1)=sigma_w_hat_adjusted(k,1)./sigma_w_hat(k,1);
%         adjustment_factor2(k,1)=sigma_p_hat_adjusted(k,1)./sigma_p_hat(k,1);

        if adjustment_factor(k,1)<0.5 | adjustment_factor(k,1) >2
            disp(['***** Adjustment factor sigma = ' num2str(adjustment_factor(k,1),'%0.2f') ' for mode ' num2str(k)]);
        end
        
    end
end

if scaletopsd
sigma_p_hat=sigma_p_hat_adjusted;
sigma_w_hat=sigma_w_hat_adjusted;
end


[Fc_out,Lc_out,Hc_out,~,S_par_out,R_par_out]=ssmod_maternglobal(lambda_hat,sigma_w_hat,p_mat*ones(size(lambda_hat)),omega_axis,tau_axis,true);










if strcmpi(doPlot,'yes')
    
plotcf(tau_axis,R_target,tau_axis,R_par_out,'type','auto','xlim',[0 max(tau_axis)],'legend',{'Target' 'Matern fit' });

plotpsd(omega_axis,S_target,omega_axis,S_par_out*2,'xlim',[0 omega_axis(end)],'type','auto','log','yes','legend',{'Target' 'Matern fit'});

end

end

