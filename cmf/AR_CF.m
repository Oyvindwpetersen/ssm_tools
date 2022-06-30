function [R_AR]=AR_CF(A_r,B,j_lag,l_lag,q)

%% Covariance function for AR-model
% From Gallego-Castillo 2021, A tutorial on reproducing a predefined autocovariance function
% through AR models: application to stationary homogeneous isotropic turbulence

% Eq 46:
% z[n] = sum_r_1toP( A[j(r)]*z[n-j(r)] ) + B*epsilon[n]
% Cov(epsilon)=E[epsilon*epsilon^T]=I

% Inputs:
% A_r: AR-matrix for selcted lags
% Sigma: input matrix
% j_lag: lags for the AR-matrices
% l_lag: desired lags for the output CF
% q: cutoff where the CF is approx zero, see discussion in Gallego Eq. 58

% Outputs:
% R_AR: 

%%

Phi_r=AR_coeffconv(A_r,'cell');

n=size(Phi_r{1},1);

% Populate full Phi, zero to zero where coefficient is not considered
Phi=cell(1,max(j_lag));
for j=1:max(j_lag)
    
    ind=find(j==j_lag);
    
    if ~isempty(ind)
        Phi{j}=Phi_r{ind};
    else
        Phi{j}=sparse(n,n);
    end

end

p=max(j_lag);

% Eq 55
Psi_0=speye(n);

Psi=cell(1,q);

for i=1:q

    Sum_temp=0;
    for j=1:min(i,p)
        
        if (i-j)>0
            Psi_i_minus_j=Psi{i-j};
        elseif (i-j)==0
            Psi_i_minus_j=Psi_0;
        end
        
        Sum_temp=Sum_temp+Phi{j}*Psi_i_minus_j;
        
    end
    
    Psi{i}=Sum_temp;    
    
end



%% Eq 58

% % Precalculate product
% BBt=B*B.';
% 
% if any(l_lag)<0
%     error('Only positive lags supported');
% end
% 
% % Precalculate product Psi_i*B*B^T
% for ind=1:q        
%    Precalc_term{ind}=Psi{ind}*BBt;
% end
% 
% R_AR=zeros(n,n,length(l_lag));
% for ind=1:length(l_lag)
%     
%     l=l_lag(ind);
%     
%     % if l>q then CF is zero by definition (approximately)
%     if l>q
%         continue
%     elseif l==0
%         
%         % if l=0 then CF is sum of Psi_i*B*B^T*Psi_i^T
%         
%         % Start sum with value at l=0 to avoid index zero, loop from l+1
%         Sum_temp=BBt;
%         for i=(l+1):q
%             Sum_temp=Sum_temp+Psi{i}*BBt*Psi{i}.';
%         end           
%         R_AR(:,:,ind)=Sum_temp;
%         continue
%     end
%     
%     Sum_temp=0;
%     for i=l:q
%         
%         if (i-l)>0
%             Psi_i_minus_l=Psi{i-l};
%         elseif (i-l)==0
%             Psi_i_minus_l=eye(n);
%         end
%         
%         Sum_temp=Sum_temp+Precalc_term{i}*Psi_i_minus_l.';
% 
%     end
%     
%     R_AR(:,:,ind)=Sum_temp;
%     
% end

%% Eq 58

% Precalculate product

if any(l_lag)<0
    error('Only positive lags supported');
end

R_AR=zeros(n,n,length(l_lag));

% Precalculate
Psi_i_B=cell(q,1);
for i=1:q
    Psi_i_B{i}=Psi{i}*B;
end

for ind=1:length(l_lag)
    
    l=l_lag(ind);
    
    % if l>q then CF is zero by definition (approximately)
    if l>q
        continue
    elseif l==0
        
        % if l=0 then CF is sum of Psi_i*B*B^T*Psi_i^T
        
        % Start sum with value at i=0 to avoid index zero, loop from l+1
        Sum_temp=B*B.';
        for i=(l+1):q
            Sum_temp=Sum_temp+Psi_i_B{i}*Psi_i_B{i}.';
        end
        
        R_AR(:,:,ind)=Sum_temp;
        continue
    end
    
    Sum_temp=0;
    for i=l:q
        
        if (i-l)>0
            Psi_i_minus_l_B=Psi_i_B{i-l};
        elseif (i-l)==0
            Psi_i_minus_l_B=speye(n)*B;
        end
        
        Sum_temp=Sum_temp+Psi_i_B{i}*Psi_i_minus_l_B.';
    end
    
	R_AR(:,:,ind)=Sum_temp;
        
end

%%

% if 0 
%     
% % Eq 58
% for ind=1:length(l_lag)
%     
%     l=l_lag(ind);
%     
%     if l>q
%         continue
%     end
%     
%     Sum_temp=0;
%     for i=l:q
%         
%         if (i-l)>0
%             Psi_i_minus_l=Psi{i-l};
%         elseif (i-l)==0
%             Psi_i_minus_l=eye(n);
%         end
%         
%         if i==0
%             Psi_i=Psi_0;
%         elseif i>0
%             Psi_i=Psi{i};
%         end
%         
%         Sum_temp=Sum_temp+Psi_i*BBt*Psi_i_minus_l.';
% 
%     end
%     
%     R_AR(:,:,ind)=Sum_temp;
%     
% end
% 
% end


