function [A_ss,B_ss]=ssmod_ar(A_block,B,j_lag)


%% PSD of AR-model

% AR-form:
% y[n] = sum_r_1toP( A_r[j(r)]*y[n-j(r)] ) + B*eta[n]

% State-space form:
% x[n]=A_ss*x[n-1]+B_ss*eta[n]

% |y[n]  |     | A1 A2 A3 A4 | |y[n-1]|   |B| eta[n]
% |y[n-1]|   = | I  0  0  0  | |y[n-2]| + |0|
% |y[n-2]|     | 0  I  0  0  | |y[n-3]|   |0|
% |y[n-3]|     | 0  0  I  0  | |y[n-4]|   |0|

% Inputs:
% A_block: cell of square A-matrices: A_block={A[j(1)],A[j(2)],...} or in 3d: A_block(:,:,1)=A[j(1)], A(:,:,2)=A[j(2)]
% B: input matrix
% j_lag: vector with lags for the AR-model. The non-sparse AR has j_lag=[0,1,2,3,...].

% Inputs:
% A_ss: state matrix
% B_ss: input matrix

%%

j_max=max(j_lag);

if iscell(A_block)
	ns=size(A_block{1},1);
	for k=1:length(A_block)
		A_mat3d(:,:,k)=A_block{k};
	end
else
	A_mat3d=A_block;
	ns=size(A_mat3d,1)
end

% Fill A_all with lags from [0,1,2,3,...,j_max]
% A_stack is zero for lags not in j_lag

A_stack=sparse(ns,ns*j_max);

for k=1:j_max
    ind=find(k==j_lag);
    if ~isempty(ind)
		ind_col=[1:ns]+(k-1)*ns;
		A_stack(1:ns,ind_col)=A_mat3d(:,:,ind);
    end
end

A_ss=[A_stack ; speye(ns*(j_max-1)) sparse(ns*(j_max-1),ns)];

B_ss=[B;sparse(size(A_ss,1)-ns,size(B,2))];
B_ss=sparse(B_ss);
