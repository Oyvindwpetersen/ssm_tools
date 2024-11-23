function [negLL,negLL_1,negLL_2,negLL_1h,negLL_2h,Nred]=loglik_evidence2(S,y,varargin)

%% Log likelihood of model evidence in state-space function
%
% Inputs:
% S: covariance matrix
% y: data time series
%
% Inputs:
% negLL: (negative) log likelihood
% negLL_1: contribution from complexity
% negLL_2: contribution from data fit
% negLL_1h: history of complexity
% negLL_2h: history of data fit
% Nred: number of data point, if cut start and end
%
%%

p=inputParser;

addParameter(p,'cut',[0 0],@isnumeric)
addParameter(p,'rejectcov',true,@islogical)
parse(p,varargin{:});

cut=p.Results.cut;
rejectcov=p.Results.rejectcov;

%%
if cut(1)>0 | cut(2)>0
    
    ind_start=[1+cut(1)];
    ind_end=size(y,2)-cut(2);
    
    range_keep=[ind_start:ind_end];
    y=y(:,range_keep);
    
    Nred=size(y,2);
    
else
    Nred=size(y,2);
end

%% Model complexity

N=size(y,2);

det_S=det(S);

% If determinant inf, then set to nan
if isinf(abs(det_S))
    det_S=NaN;
    
    negLL_1=NaN; negLL_2=NaN; negLL_2=NaN; negLL_2h=NaN; negLL=NaN;
    return
    
end
    
negLL_1=N*logdet(Omega);

negLL_1h=negLL_1/N;

%% Data fit

if any(any(isnan(S)))
    S
    error('S contains nan');
end

if any(any(isinf(S)))
    S
    error('S contains inf');
end

Omega=forcesym(S);
[V,D]=eig(Omega);

Omega_half=V*D.^0.5;
Omega_half_inv=eye(size(Omega_half))/Omega_half;
Temp2=Omega_half_inv*y;
negLL_2_temp=sum(Temp2.^2,1);
negLL_2=sum(negLL_2_temp);

negLL_2h=negLL_2_temp;

if any(abs(diag(D))>1e100)
    
    negLL_1=NaN; negLL_2=NaN; negLL_2=NaN; negLL_2h=NaN; negLL=NaN;
    return
    
end

if any(diag(D)<-1e-6)
    
    negLL_1=NaN; negLL_2=NaN; negLL_2=NaN; negLL_2h=NaN; negLL=NaN;
    return
    
end

% negLL_2_check=0;
% for k=1:length(e)
%     negLL_2_check=negLL_2_check+e(:,k).'/Omega*e(:,k);
% end

%% Sum

negLL=negLL_1+negLL_2;

