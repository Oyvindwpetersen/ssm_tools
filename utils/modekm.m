function modekm(K,M)

%%

K=mod.Ke
M=mod.Me


[v,d]=eig(K,M);

omega=sqrt(diag(d));

alpha=1./sqrt(diag(v.'*M*v));

phi=NaN*ones(size(v));

for k=1:size(K,1)
    phi(:,k)=v(:,k)*alpha(k);
end
