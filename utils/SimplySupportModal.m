function [Kg,Mg,f,phi,phi_full,K,M,K_red,M_red,S_red]=SimplySupportModal(L,E,I,m,N_el)

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
    
% N_el=100
N_node=N_el+1;

K=zeros(N_node*2,N_node*2);
M=zeros(N_node*2,N_node*2);

for k=1:N_el
    
    L_el=L/N_el;
    ind=[1:4]+(k-1)*2;
    K(ind,ind)=K(ind,ind)+K_mat(E,I,L_el);
    M(ind,ind)=M(ind,ind)+M_mat(m,L_el);
    
end

%% Eliminate BCs

ind_keep=setdiff(1:length(K),[1 length(K)-1]);
K_red=K(ind_keep,ind_keep);
M_red=M(ind_keep,ind_keep);

%% Modes

[v,d]=eig(K_red,M_red);

dd=diag(d);
[~,i_sort]=sort(dd);
d_sort=dd(i_sort);
v=v(:,i_sort);

% w=sqrt(dd(i_sort))
% f=w/(2*pi)

% f_ref=1/(2*pi)*([1:10]*pi/L).^2*sqrt(E*I/rho)

% figure(); hold on;
% plot(f(1:10)./f_ref.','-ob');

%% Scale modes

omega=sqrt(d_sort);
f=omega./(2*pi);

Mg_temp=v.'*M_red*v;

factor(1,:)=diag(Mg_temp);

phi=v./sqrt(factor);

Kg=phi.'*K_red*phi;
Mg=phi.'*M_red*phi;

Kg=sparse(diag(diag(Kg)));
Mg=sparse(diag(diag(Mg)));

S_red=sparse(ind_keep,1:length(ind_keep),ones(1,length(ind_keep)),length(K),length(ind_keep));

phi_full=S_red*phi;
