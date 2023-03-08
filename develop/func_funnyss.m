function f5=func_funnyss(a,k0,m1,m2,m3)


%% 

a=1;
k0=40;
k1=40;
m1=0;
m2=0
m3=00;

%
L=30;
E=210e9;
I=1e-8;
m=50;
N_el=100;

[Kg,Mg,f,phi,phi_full,K0,M0,K_red,M_red,S_red]=SimplySupportModal(L,E,I,m,N_el);

K=K0;
M=M0;

p=polyfit([-1 0 1],[k0 k1 k0],2);

k_spring=polyval(p,[-1 0 1]);

x=linspace(-1,1,101);

k_spring=polyval(p,x);

% k_spring=(1*x.^2+a)*k0;

for k=1:1:101
    i=k*2-1;
    K(i,i)=K(i,i)+k_spring(k);
end

x_node=linspace(0,1,101);

m_spring=[m1 m2 m3 m2 m1];

no=0;
for k=1:20:100
    i=k*2-1;
    
    no=no+1;
    M(i,i)=M(i,i)+m_spring(no);
    
end

ind_keep=1:length(K);
% ind_keep=setdiff(1:length(K),[1 length(K)-1]);
K_red=K(ind_keep,ind_keep);
M_red=M(ind_keep,ind_keep);

[v,d]=eig(K_red,M_red);

dd=diag(d);
[~,i_sort]=sort(dd);
d_sort=dd(i_sort);
v=v(:,i_sort);

f=sqrt(d_sort);

% v_full=S_red*v;
v_full=v;


f5=f(1:10);

f5=f5./f5(1);

f5

%%

close all

figure(); hold on; grid on;
plot(f(1:20),'o');

figure(); hold on; grid on;

for k=1:5
    plot(v_full(1:2:end,k));
end

tilefigs([2 2],'l');