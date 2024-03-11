function [lambda_export,phi_export,w_export,xi_export,psi_export]=eigAC(A,C,dt,uniqueOnly)

%%

if ~iscell(A)
    A={A};
    C={C};
   exportCell=false; 
else
   exportCell=true; 
end

if nargin==3
    uniqueOnly='';
end

%%

for j=1:length(A)

[v,d_discrete]=eig(A{j});
d_cont=log(d_discrete)./dt;

[~,i_sort]=sort(abs(diag(d_cont)));
d_cont=diag(d_cont); d_cont=d_cont(i_sort);
v=v(:,i_sort);

v=v./max(abs(v),[],1);

lambda=d_cont;
phi=v;
w=abs(lambda);
xi=-real(lambda)./abs(lambda);

psi=C{j}*phi;

for k=1:size(psi,2)
   [~,rotrad,~,~]=mpcweighted(psi(:,k));
	psi(:,k)=conj(psi(:,k))*exp(1i*rotrad);
	phi(:,k)=conj(phi(:,k))*exp(1i*rotrad);
end

if strcmpi(uniqueOnly,'unique');
    [w,indExport]= uniquetol(w,1e-6);
else
    indExport=1:length(lambda);
end    
    
lambda_export{j}=lambda(indExport);
phi_export{j}=phi(:,indExport);
w_export{j}=w;
xi_export{j}=xi(indExport);
psi_export{j}=psi(:,indExport);

end

%%

if exportCell==false
    lambda_export=lambda_export{1};
    phi_export=phi_export{1};
    w_export=w_export{1};
    xi_export=xi_export{1};
    psi_export=psi_export{1};
end


