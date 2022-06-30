function checkDim(A,rows,cols);

%%

[n1,n2]=size(A);

if n1~=rows & n2==cols
    error(['Row dimension of matrix wrong']);
elseif n1==rows & n2~=cols
    error(['Column dimension of matrix wrong']);   
elseif n1~=rows & n2~=cols
    error(['Row and column dimension of matrix wrong']);
end
    
