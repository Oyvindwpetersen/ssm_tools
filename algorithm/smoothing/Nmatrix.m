function N_L=Nmatrix(A,G,L)
%%

ns=size(A,1);
nd=size(G,1);

N0=zeros(nd,ns);
N1=[zeros(nd,ns) ; G];

N_k=cell(1,L);
N_k{1}=N1;

for k=2:L
    
    Obs_kminus=ExtendedObservability(A,G,k-1);
    N_k{k}=zeros((k+1)*nd,k*ns); 
    
    [i1,i2]=size(N_k{k-1});
    [j1,j2]=size(Obs_kminus);    
    
    N_k{k}((end-i1+1):end,:)=[Obs_kminus N_k{k-1}];
end

if L==0
    N_L=N0;
elseif L==1
    N_L=N1;
elseif L>=2
    N_L=N_k{L};
end
    
