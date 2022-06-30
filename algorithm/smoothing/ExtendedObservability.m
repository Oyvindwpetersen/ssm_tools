function Obs_L=extendedObservability(A,G,L);
%%

ns=size(A,1);
nd=size(G,1);

Obs_L=zeros((L+1)*nd,ns);

index_range=[1:nd];
for r=1:(L+1)
	Obs_L(index_range,:)=G*A^(r-1);
	index_range=index_range+nd;
end

