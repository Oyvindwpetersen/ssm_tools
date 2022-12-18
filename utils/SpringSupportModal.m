function [Kg,Mg,f,phi,phi_full,K,M,K_red,M_red,S_red]=SpringSupportModal(L,E,I,m,N_el,N_modes,k0,k_end)

%% 

% k0=40;
%
% L=30;
% E=210e9;
% I=1e-8;
% m=50;
% N_el=100;

[~,~,~,~,~,K0,M0,~,~,~]=SimplySupportModal(L,E,I,m,N_el);

K=K0;
M=M0;

k1=k0;
p=polyfit([-1 0 1],[k0 k1 k0],2);

k_spring=polyval(p,[-1 0 1]);

N_node=N_el+1;
x=linspace(-1,1,N_node);
k_spring=polyval(p,x);

for k=1:1:(N_el+1)
    i=k*2-1;
    K(i,i)=K(i,i)+k_spring(k);
end

K(1,1)=K(1,1)+k_end;
K(end-1,end-1)=K(end-1,end-1)+k_end;

ind_keep=1:length(K);
% ind_keep=setdiff(1:length(K),[1 length(K)-1]);

K_red=K(ind_keep,ind_keep);
M_red=M(ind_keep,ind_keep);

[v,d]=eig(K_red,M_red);

dd=diag(d);
[~,i_sort]=sort(dd);
d_sort=dd(i_sort);
v=v(:,i_sort);

omega=sqrt(d_sort);

% v_full=S_red*v;
% v_full=v;

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

phi_full=phi;
