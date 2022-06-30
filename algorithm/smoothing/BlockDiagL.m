function [X_L_i]=BlockDiagL(X,L,i)

%% 

[n1,n2]=size(X);

% Number of matrices
n_matrices=L-abs(i);

% Fill the shifted diagonal (i'th below main diagonal) with ones
block_index=diag(ones(1,n_matrices),-i);
block_index(block_index==0)=2;

% Find the index of positions to fill
[i1,i2]=find(block_index==1);

% Fill the big matrix on the shifted diagonal
X_L_i=sparse(n1*L,n2*L);
for k=1:n_matrices
    range1=(i1(k)-1)*n1+[1:n1];
    range2=(i2(k)-1)*n2+[1:n2];
    
    X_L_i(range1,range2)=X;
end

