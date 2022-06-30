function [Q_L,R_Lplus,S_L]=extendedNoiseCovariance(Q,R,S,L,i)

%%

ns=size(Q,1);
nd=size(R,1);

Q_L=BlockDiagL(Q,L,i);

R_Lplus=BlockDiagL(R,L+1,i);

if i>=0
    S_temp=BlockDiagL(S,L,i);
    S_L=[S_temp zeros(L*ns,nd)];
elseif i<0
	S_temp=BlockDiagL(S,L,i+1);
    S_L=[zeros(L*ns,nd) S_temp];
end  

if util(Q_L)<1e-4
    Q_L=sparse(Q_L);
end

if util(R_Lplus)<1e-4
    R_Lplus=sparse(R_Lplus);
end

if util(S_L)<1e-4
    S_L=sparse(S_L);
end
        

checkDim(Q_L,L*ns,L*ns);

checkDim(R_Lplus,(L+1)*nd,(L+1)*nd);

checkDim(S_L,L*ns,(L+1)*nd);


