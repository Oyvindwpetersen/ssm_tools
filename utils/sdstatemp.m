function sd_emp=sdstatemp(x_ref,x_test,n_cut)

%% Calculate statistics (SD) of time series estimate(s) vs. reference
%
% Inputs:
% x_ref: matrix with each signal as row vector
% x_test: matrix (or cell if multiple tests) with each signal as row vector
% n_cut: [1,2] matrix with number of samples to trim at start and end 
%
% Outputs:
% plot_std_emp: cell with SDs 
%
%%

% Ensure cell
if ~iscell(x_test)
    x_test={x_test};
end

% Default no trimming of 
if nargin<3
    n_cut=[];
end

if isempty(n_cut)
    n_cut=[0 0];
end

% If one number given, apply to start and end
if length(n_cut)==1
    n_cut=[n_cut n_cut];
end

t_range=[(1+n_cut(1)):(size(x_ref,2)-n_cut(2))];

% Calculate std
sd_emp={};
for k=1:length(x_test)

    sd_emp{k}=[std(x_ref(:,t_range)-x_test{k}(:,t_range),0,2) ];

end

