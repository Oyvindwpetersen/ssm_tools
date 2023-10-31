function [Kg,Mg,f,phi,K,M,phi_red,K_red,M_red,S_red,S_red_inv,nodecoord,doflabel,doflabel_red]=simplysupportbeam(L,E,I,m,N_el,varargin)

%%

% Inputs: 
% 
% 

% Outputs: 
% 

%% Input handling

p=inputParser;
p.KeepUnmatched=true;

addParameter(p,'bc','support',@ischar) %springsupport
addParameter(p,'k_end',0,@isnumeric)
parse(p,varargin{:})

bc=p.Results.bc;
k_end=p.Results.k_end;

%%
N_node=N_el+1;

nodelabel=1:N_node;
xcoord=linspace(0,L,length(nodelabel)); ycoord=xcoord*0; zcoord=xcoord*0;
nodecoord=[nodelabel.' xcoord.' ycoord.' zcoord.'];

doflabel=genlabels({'U' 'UR'},nodelabel);

%% Element matrices

K_mat=@(E,I,L) E*I/L^3*...
    [ 12 6*L -12 6*L ;
      6*L 4*L^2 -6*L 2*L^2 ; 
      -12 -6*L 12 -6*L ; 
      6*L 2*L^2 -6*L 4*L^2];

M_mat=@(m,L) m*L/420*...
    [156 22*L 54 -13*L ;
     22*L 4*L^2 13*L -3*L^2 ;
     54 13*L 156 -22*L ;
     -13*L -3*L^2 -22*L 4*L^2];
	 
%% Build global matrices

K=zeros(N_node*2,N_node*2);
M=zeros(N_node*2,N_node*2);

for k=1:N_el
    
    L_el=L/N_el;
    ind=[1:4]+(k-1)*2;
    K(ind,ind)=K(ind,ind)+K_mat(E,I,L_el);
    M(ind,ind)=M(ind,ind)+M_mat(m,L_el);
    
end

%% Assign BCs

if strcmpi(bc,'support')

    ind_keep=setdiff(1:length(K),[1 length(K)-1]);
    S_red=sparse(ind_keep,1:length(ind_keep),ones(1,length(ind_keep)),length(K),length(ind_keep));
    S_red_inv=pinv(full(S_red));
    
    K_red=K(ind_keep,ind_keep);
    M_red=M(ind_keep,ind_keep);
    
elseif strcmpi(bc,'spring')
    
    K(1,1)=K(1,1)+k_end;
    K(end-1,end-1)=K(end-1,end-1)+k_end;
    
    S_red=eye(size(K));
    S_red_inv=S_red;
    
    K_red=K;
    M_red=M;
    
end

%% Modes

K_red=S_red.'*K*S_red;
M_red=S_red.'*M*S_red;

[v,d]=eig(K_red,M_red);

dd=diag(d);
[~,i_sort]=sort(dd);
d_sort=dd(i_sort);
v=v(:,i_sort);

%% Scale modes

omega=sqrt(d_sort);
f=omega./(2*pi);

Mg_temp=v.'*M_red*v;

factor(1,:)=diag(Mg_temp);

phi_red=v./sqrt(factor);

Kg=phi_red.'*K_red*phi_red;
Mg=phi_red.'*M_red*phi_red;

Kg=sparse(diag(diag(Kg)));
Mg=sparse(diag(diag(Mg)));

phi=S_red*phi_red;

% w=sqrt(dd(i_sort))
% f=w/(2*pi)
