%% Time idare time

% S vs no S
% scaling vs noscaling

% Conc:
% S or no S plays no role
% noscaling gives 1/3 reduction in computational time

clc

F=eye(40)+randn(40)/10000;

H=randn(20,40);

Q=eye(40)*100;
R=eye(20)+randn(20); R=R*R.';
S=zeros(40,20);

S=randn(size(S))/1000;

tic
for k=1:100
[P_k_kmin_ss1,~,~,info]=idare((F).',(H).',Q,R,S,'noscaling');
end
toc

tic
for k=1:100
[P_k_kmin_ss2,~,~,info]=idare((F).',(H).',Q,R,S);
end
toc

tic
for k=1:100
[P_k_kmin_ss3,~,~,info]=idare((F).',(H).',Q,R,'noscaling');
end
toc

tic
for k=1:100
[P_k_kmin_ss4,~,~,info]=idare((F).',(H).',Q,R);
end
toc