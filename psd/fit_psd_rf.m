function [n_opt,d_opt,alpha_opt,S_opt,rn_opt,rd_opt]=fit_psd_rf(omega,S_target,order_n,order_d,varargin)

%% Fit rational function to target spectral density
%
% S(omega)=alpha*N(omega)/D(omega)
%
% N(omega)=n0+n2*omega^2+n4*omega^4+...+n_2M*omega^(2M)
% D(omega)=d0+d2*omega^2+d4*omega^4+...+d_2K*omega^(2K)
%
% The highest order coefficients (n_2M and d_2K) are set to 1 by definition
%
% Inputs:
% omega: frequency vector
% S_target: target spectral density (two-sided)
% order_n: order of N-polynomial, must be even
% order_d: order of D-polynomial, must be even
% varargin: see below
%

%% Inputs

p=inputParser;
p.KeepUnmatched=true;
addParameter(p,'plot',true,@islogical) % plot spectrum
addParameter(p,'forcezero',false,@islogical) % force to origin, S(omega=0)=0, equivalent to n0=0
addParameter(p,'alpha0',[],@isnumeric) % initial guess for alpha
addParameter(p,'n_ini',[],@isnumeric) % initial guess for coefficients n=[n_2M ... n2 n0]
addParameter(p,'d_ini',[],@isnumeric) % initial guess for coefficients d=[d_2K ... d2 d0]
addParameter(p,'rn_ini',[],@isnumeric) % initial guess for roots of N(omega)
addParameter(p,'rd_ini',[],@isnumeric) % initial guess for  roots of D(omega)
addParameter(p,'penaltylog',false,@islogical)

parse(p,varargin{:})
do_plot=p.Results.plot;
forcezero=p.Results.forcezero;
alpha0=p.Results.alpha0;
n_ini=p.Results.n_ini;
d_ini=p.Results.d_ini;
rn_ini=p.Results.rn_ini;
rd_ini=p.Results.rd_ini;
penaltylog=p.Results.penaltylog;

%% Input check

[S_target]=inputcheck(omega,S_target,order_n,order_d,n_ini,d_ini,forcezero);

%% 

% Number of real and imag
no_z_im=ceil(order_n/2/2);
no_z_re=order_n/2-no_z_im;

% If N(omega=0)=0, then two roots must be set explicitly to zero
% This is obtained by removing two roots, then adding [0 ; 0] later
if forcezero  & any(order_n==[4 8 12 16 20 24 28])
    no_z_re=no_z_re-1;
end

if forcezero  & any(order_n==[2 6 10 14 18 22 26])
    no_z_im=no_z_im-1;
end

if any(order_n>[28])
    error('Such high orders not implemented, code above this line must be fixed');
end

no_p_im=ceil(order_d/2/2);
no_p_re=order_d/2-no_p_im;

% Determine ranges for each parameter in x-vector, x=[alpha ; z real ; z imag ; p real ; p imag];
offset=1;
range_n_re=1+[1:no_z_re];

if ~isempty(range_n_re)
	offset=range_n_re(end);
end
range_n_im=offset+[1:no_z_im];

if ~isempty(range_n_im)
    offset=range_n_im(end);
end
range_a_re=offset(end)+[1:no_p_re];

if ~isempty(range_a_re)
    offset=range_a_re(end);
end
range_a_im=offset+[1:no_p_im];

% Function producing spectral density from roots
psd_rf_func=@(x) psd_rf_root(omega,x(1),x(range_n_re),x(range_n_im),x(range_a_re),x(range_a_im),forcezero);

%% Default initial values for poles

z_re0=linspace(0.25,1,no_z_re).';
z_im0=linspace(0.25,0.25,no_z_im).'/10;

[~,ind_max]=max(S_target);

spread=[0 0.1 0.2 0.3 0.4 0.5];
p_re0=linspace((1-spread(no_p_re))*omega(ind_max),(1+spread(no_p_re))*omega(ind_max),no_p_re).';

p_im0=linspace(0.5,0.5,no_p_im).';

%% Override initial values, if provided

% Initial values for polynomial coeffcients, even terms only
if ~isempty(n_ini) & ~isempty(d_ini)

    [rn_ini,z_re0,z_im0]=rooteven(n_ini);
    [rd_ini,p_re0,p_im0]=rooteven(d_ini);

% Initial values for roots, must be conjugate pairs
elseif ~isempty(rn_ini) & ~isempty(rd_ini)
    
    n_ini_all=poly(rn_ini); n_ini=n_ini_all(1:2:end); n_ini_odd=n_ini_all(2:2:end);
    d_ini_all=poly(rd_ini); d_ini=d_ini_all(1:2:end); d_ini_odd=n_ini_all(2:2:end);

    if any(abs(n_ini_odd)>1e-6)
        n_ini
        n_ini_odd
        error('Odd coefficients should be zero, this indicates that the roots are not complex conjugate pairs');
    end

    if any(abs(d_ini_odd)>1e-6)
        d_ini
        d_ini_odd
        error('Odd coefficients should be zero, this indicates that the roots are not complex conjugate pairs');
    end
    
    
    [rn_ini,z_re0,z_im0]=rooteven(n_ini);
    [rd_ini,p_re0,p_im0]=rooteven(d_ini);
    
end

%% All initial values

x0=[1 ; z_re0 ; z_im0 ; p_re0 ; p_im0];

%% Alpha scaling

S0_temp=psd_rf_func(x0);

% If initial value for alpha is not provided, scale so that max(S0)=max(S_target)
if isempty(alpha0)
    x0(1)=max(S_target)./max(S0_temp)*x0(1);
else
    x0(1)=alpha0;
end

S0=psd_rf_func(x0);

% plotpsd(omega,S_target,S0,'xlim',[0 10]);

%% Scaling of objective function for easier interpretation

domega=diff(omega(1:2));

N_omega=length(omega);
target_var=sqrt( 1/(domega*N_omega)*sum(S_target.^2));

if max(abs(diff(omega)-domega))./domega > 1e-6
    warning('omega vector not evenly spaced, setting in scaling equal to 1 in least squares');
	N_omega=1;
    domega=1;
    target_var=1;
end

%% Objective function 

% Bounds
lb=0*ones(size(x0)); lb(1)=max(max(S_target)*1e-4,eps*10);
ub=20*ones(size(x0)); ub(1)=max(S_target)*100;

% Difference
func_delta= @(x) S_target-psd_rf_func(x);

% If penalized on log basis
if penaltylog
    func_delta= @(x) log10(S_target)-log10(psd_rf_func(x));
end

% Least squares objective function, normalized by variance
func_obj_nrms= @(x) sqrt( 1/(N_omega*domega)*sum(func_delta(x).^2) ) /target_var*100;

% Shift violations of bounds slightly into domain
ind_vio_lb=~(x0>lb); ind_vio_ub=~(x0<ub);
x0(ind_vio_lb)=lb(ind_vio_lb)+1e-6*(ub(ind_vio_lb)-lb(ind_vio_lb));
x0(ind_vio_ub)=ub(ind_vio_ub)-1e-6*(ub(ind_vio_ub)-lb(ind_vio_ub));

F0=func_obj_nrms(x0);

%% Initial DLAS

% PerturbSigma0=0.05*(ub-lb);
% PerturbSigma0(1)=PerturbSigma0(1)/10;
% 
% options=DLASopt('ShowText',true,'IntervalText',100,'x_lb',lb,'x_ub',ub,'L',10,'S0',[x0],...
%                 'PerturbType','randn','halflife',[],'IterationMax',20e3,...
%                 'PerturbSigma0',PerturbSigma0);
% 
% [x_star,F_star]=DLAS(func_obj_nrms,options);

%% fmincon

Aeq=[]; beq=[];

options = optimoptions(@fmincon,'Display','iter','StepTolerance',1e-6,'MaxFunctionEvaluations',2000,'ScaleProblem',false);
[x_opt1,f_opt1]=fmincon(func_obj_nrms,x0,[],[],Aeq,beq,lb,ub,[],options);

alpha_opt1=x_opt1(1);
[S_opt1,rn_opt1,n_opt1,rd_opt1,d_opt1]=psd_rf_root(omega,x_opt1(1),x_opt1(range_n_re),x_opt1(range_n_im),x_opt1(range_a_re),x_opt1(range_a_im),forcezero);

% S_opt1=psd_rf_func(x_opt1);

%% Global search

warning('off','MATLAB:singularMatrix');

problem = createOptimProblem('fmincon','x0',x_opt1,'objective',func_obj_nrms,'lb',lb,'ub',ub);
gs = GlobalSearch('XTolerance',1e-6,'StartPointsToRun','bounds','NumTrialPoints',300,'NumStageOnePoints',100,'Display','iter');
[x_opt2,f_opt2] = run(gs,problem);

alpha_opt2=x_opt2(1);
[S_opt2,rn_opt2,n_opt2,rd_opt2,d_opt2]=psd_rf_root(omega,x_opt2(1),x_opt2(range_n_re),x_opt2(range_n_im),x_opt2(range_a_re),x_opt2(range_a_im),forcezero);

warning('on','MATLAB:singularMatrix');

%%

if do_plot
    
    plotopt=struct();
    plotopt.xlim=[0 10];
    plotopt.xlabel={'Frequency [rad/s]'};
    plotopt.displayname={'Target' 'Initial' 'Optimized'};
    plotopt.lineStyle={'-' '-' '--'};
    
    plotpsd(omega,S_target,S0,S_opt2,plotopt);
end

%%

n_opt=n_opt2;
d_opt=d_opt2;
alpha_opt=alpha_opt2;
S_opt=S_opt2;
rn_opt=rn_opt2;
rd_opt=rd_opt2;

end

%% Ensure proper inputs

function [S_target]=inputcheck(omega,S_target,order_n,order_d,n_ini,d_ini,forcezero)

% Ensure 3d
if size(S_target,3)==1
    S_target_3d(1,1,:)=S_target;
    S_target=S_target_3d; clear S_target_3d
end

% Check non-negative
if any(S_target<0)
    
    ind_neg=find(S_target<0);
    omega_neg=omega(ind_neg)
    S_target_neg=S_target(ind_neg)
    
    error('Target spectral density must be non-negative');
end 

% Check inf
if any(isinf(S_target))
    
    ind_inf=find(isinf(S_target));
    omega_inf=omega(ind_inf)
    S_target_inf=S_target(ind_inf)
    
    error('Target spectral density must be non-infite');
end

% Check real

ratio=sqrt(sum(imag(S_target).^2)./sum(real(S_target).^2));

if ratio>1e-6  
    error('Target spectral density must be non-imaginary');
else
    S_target=real(S_target);
end 

if mod(order_n,2)==1
    order_n
    error('Order of nominator must be an even number');
end

if mod(order_d,2)==1
    order_d
    error('Order of denominator must be an even number');
end

if order_n>=order_d
    order_n
    order_d
    error('Order of nominator must be smaller than order of denominator');
end

if order_n==0 & forcezero
    order_n
    error('Cant have order_n=0 when forcezero is true, this would give N(omega)=n0=0, so S(omega)=0');
end

if order_d==0
    order_d
    error('Cant be zero');
end

if ~isempty(n_ini) & length(n_ini)~=(order_n/2+1)
    order_n
    n_ini
    error('Length of n_ini does not match order_n');
end

if ~isempty(d_ini) & length(d_ini)~=(order_d/2+1)
    order_d
    d_ini
    error('Length of d_ini does not match order_d');
end


if ~isempty(n_ini) & abs(n_ini(1)-1)>1e-6
    n_ini
    error('First coefficient in n_ini must be 1');
end

if ~isempty(d_ini) & abs(d_ini(1)-1)>1e-6
    d_ini
    error('First coefficient in d_ini must be 1');
end

end

