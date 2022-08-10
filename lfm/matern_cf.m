function R_out=matern_cf(tau_axis,lambda_vec,sigma_w_vec,p_vec,output_type)

%% Matern model covariance function
%
% Inputs:
% tau_axis: time axis in [s]
% lambda_vec: vector with hyperparameters, inverse length scale
% sigma_w_vec: vector with hyperparameters, standard deviations
% p_vec: vector with order
% output_type: '2d' or '3d'
%
% Outputs:
% R_out: matrix with CF
%
%%

if any(tau_axis<0)
    error('Only positive time lags allowed');
end

if nargin<3
    output_type='3d';
end

R_par_2d=zeros(length(p_vec),length(tau_axis));

for j=1:length(p_vec)
    
        if p_vec(j)>3
			error('Function only supported for p=0,1,2,3');
		end
    
        if p_vec(j)==0
            R_par_2d_temp=sigma_w_vec(j).^2 ./ (2*lambda_vec(j)).*exp(-lambda_vec(j).*tau_axis);
        elseif p_vec(j)==1
            R_par_2d_temp=sigma_w_vec(j).^2 ./ (4*lambda_vec(j).^3).*exp(-lambda_vec(j).*tau_axis).*(1+lambda_vec(j).*tau_axis);
        elseif p_vec(j)==2
            R_par_2d_temp=sigma_w_vec(j).^2 ./ (16*lambda_vec(j).^5).*exp(-lambda_vec(j).*tau_axis).*(lambda_vec(j).^2.*tau_axis.^2+3*lambda_vec(j).*tau_axis+3);
        elseif p_vec(j)==3
            R_par_2d_temp=sigma_w_vec(j).^2 ./ (96*lambda_vec(j).^7).*exp(-lambda_vec(j).*tau_axis).*(lambda_vec(j).^3*tau_axis.^3+6*lambda_vec(j).^2.*tau_axis.^2+15*lambda_vec(j).*tau_axis+15);
        end
		R_par_2d(j,:)=R_par_2d_temp;
end

if strcmpi(output_type,'3d')
	R_out=diag3d(R_par_2d);
elseif strcmpi(output_type,'2d')
	R_out=R_par_2d;
end