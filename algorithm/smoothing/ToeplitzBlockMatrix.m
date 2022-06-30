function H=ToeplitzBlockMatrix(A,B,G,J,L)
%%

C=cell(1,L);

C{1}=J;
for r=2:(L+1)
    C{r}=G*A^(r-2)*B;    
end

n=length(C);

[s1,s2]=size(J);
C{n+1}=zeros(s1,s2);

toeplitz_index=toeplitz(1:n);

for i=1:(n-1)
    for j=(i+1):n
        toeplitz_index(i,j)=n+1;
    end
end

C_stack=C(toeplitz_index);
H = cell2mat(C_stack);


