function [x_filt p_filt Pkk Ppkk] = JIS(A,B,G,J,y,x0,xdot0,R,Q,P01,varargin)


p=inputParser;
addParameter(p,'inverse','conv',@ischar)
addParameter(p,'gainsupplied','no',@ischar)
addParameter(p,'gain',NaN,@isstruct)
parse(p,varargin{:});
inverse = p.Results.inverse;
gainsupplied = p.Results.gainsupplied;
gain = p.Results.gain;

%%

% A=fid.A;
% B=fid.B;
% G=fid.Gs;
% J=fid.Js;
% x0=fid.x0;
% xdot0=fid.xdot0;
% y=fid.ys;
% P01=fid.P01;
% S=fid.S;
% R=fid.Rs;
% Q=fid.Q;
% p0=fid.p0;
% Pp0=fid.Pp0;
% P0=fid.P0;

ns=size(A,1);
nt=length(y);
nd=size(G,1);
np=size(B,2);

%zero matrices
x_pred=zeros(ns,nt);
x_filt=zeros(ns,nt);
p_filt=zeros(np,nt);
innov=zeros(nd,nt);

%assign initial values
x_pred(:,1)=[x0;xdot0];
Pkk_=P01;

%
lambdaRk=zeros(nd,nt);
lambdaJRkJ=zeros(np,nt);
lambdaJPpkkJ=zeros(nd,nt);


%% Gain supplied as state state values?
if  strcmpi(gainsupplied,'yes')
    tstart=tic;
    if isnumeric(gain)
    if isnan(gain)
        disp('Error in gain: NaN')
        return
    end
    end

    Rk=gain.Rk;
    Lk=gain.Lk;
    Mk=gain.Mk;
    Ppkk=gain.Ppkk;
    Pkk=gain.Pkk;
    Pkk_=gain.Pkk_;

    for k=1:nt;

        %input estimation
        p_filt(:,k)=Mk*(y(:,k)-G*x_pred(:,k));

        %measurement update
        x_filt(:,k)=x_pred(:,k)+Lk*(y(:,k)-G*x_pred(:,k)-J*p_filt(:,k));

        %time update
        x_pred(:,k+1)=A*x_filt(:,k)+B*p_filt(:,k);

    end
    
    telapsed=toc(tstart);
    disp(['JIS excluding gain calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed*10^3./nt) ' seconds per 1k steps']);
    return
end

%% Conventional one step recursive
tstart=tic;
for k=1:nt;

    %input estimation
    if strcmp(inverse,'conv')
        Rk =G*Pkk_*G'+R;  Rk=forcesym(Rk);
        Mk = (J'/Rk*J)\J'/Rk;
        p_filt(:,k)=Mk*(y(:,k)-G*x_pred(:,k));
        Ppkk=eye(np)/(J'/Rk*J);  Ppkk=forcesym(Ppkk);

    elseif  strcmp(inverse,'svd')
        Rk =G*Pkk_*G'+R;  Rk=forcesym(Rk);
        a2=svdinv(Rk,1e-12);
        a1=svdinv(J'*a2*J,1e-12);
        Mk = a1*J'*a2;
        p_filt(:,k)=Mk*(y(:,k)-G*x_pred(:,k));
        Ppkk=a1; Ppkk=forcesym(Ppkk);
    end
    
    %measurement update
    Lk=Pkk_*G'/Rk;
    x_filt(:,k)=x_pred(:,k)+Lk*(y(:,k)-G*x_pred(:,k)-J*p_filt(:,k));
    Pkk = Pkk_ - Lk*(Rk-forcesym(J*Ppkk*J'))*Lk'; Pkk=forcesym(Pkk);
    Pxpkk=-Lk*J*Ppkk;
    Ppxkk=Pxpkk';

    %time update
    x_pred(:,k+1)=A*x_filt(:,k)+B*p_filt(:,k);
    Pkk_ = [A B]*[ Pkk Pxpkk;Ppxkk Ppkk]*[A';B']+Q; Pkk_=forcesym(Pkk_);
    
	innov(:,k)=(y(:,k)-G*x_pred(:,k)-J*p_filt(:,k));

    %eigenvalues
    
%          c1(k)=cond(Rk);
%          c2(k)=cond(J'/Rk*J);
%          c3(k)=cond(J*Ppkk*J');
         
%          
%         [V1,lambda_orig1] = eig(Rk);
%         [lambda1,I1] = sort(diag(lambda_orig1),1,'descend');
%         lambdaRk(:,k) = lambda1;
%         
%         [V2,lambda_orig2] = eig(forcesym(J'/Rk*J));
%         [lambda2,I2] = sort(diag(lambda_orig2),1,'descend');
%         lambdaJRkJ(:,k) = lambda2;
%         
%         [V3,lambda_orig3] = eig(forcesym(J*Ppkk*J'));
%         [lambda3,I3] = sort(diag(lambda_orig3),1,'descend');
%         lambdaJPpkkJ(:,k) = lambda3;
%         
end
telapsed=toc(tstart);
disp(['JIS calculated in ' sprintf('%2.1f', telapsed) ' seconds, ' sprintf('%2.1f', telapsed*10^3./nt) ' seconds per 1k steps']);



%%
return
close all

% bar3z(abs(Lk)./repmat(max(abs(Lk),[],2),1,fid.nd),'name2',fid.a_o,'name1',fid.modelabel)

figure(); hold on; grid on;
plot(lambdaRk');
set(gca,'YScale','log');

figure(); hold on; grid on;
plot(lambdaJRkJ');
set(gca,'YScale','log');

figure(); hold on; grid on;
plot(lambdaJPpkkJ');
set(gca,'YScale','log');
% 
% figure(); hold on; grid on;
% plot(abs(lambdaJPpkkJ(1,:)./lambdaJPpkkJ(end,:)));
% set(gca,'YScale','log');

%%
% return
close all

figure(); hold on; grid on;
plot(c1)
set(gca,'YScale','log');
title('Cond Rk')

figure(); hold on; grid on;
plot(c2)
set(gca,'YScale','log');
title('Cond J^T Rk J')

figure(); hold on; grid on;
plot(c3)
set(gca,'YScale','log');


