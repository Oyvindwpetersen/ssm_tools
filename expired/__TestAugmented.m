%%

addPathTemp2('C:\Cloud\OD_OWP\Work\Hardanger\GPLFM\latentforcemodel');

Aa=[mod.A mod.B ; zeros(mod.np,mod.nm*2) eye(mod.np)];

Ga=[mod.G mod.J];

Qp_val=10.^[-10:5];

clear p_all e_norm p_norm objf_1 objf_2
for k=1:length(Qp_val)

    Qp=Qp_val(k)*diag([1 1]);
    Q=blkdiag(fid.Q,Qp);
    R=fid.R;
    S=zeros(size(R,1),size(Q,1)).';
    x0=[fid.x0 ; zeros(mod.np,1)];
    P_0_0=eye(size(Aa));
    y=fid.y;
    
    [xa_est,xa_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KalmanFilter(Aa,Ga,Q,R,S,y,x0,P_0_0);
    
    e=fid.y_clean-Ga*xa_est;
    x=xa_est(1:end-mod.np,:);
    p=xa_est(end-mod.np+1:end,:);
    
    e_norm(k)=sum(sum(e.^2,1),2);
    p_norm(k)=sum(sum(p.^2,1),2);
    
    x_all{k}=x;
    p_all{k}=p;
    
    p_err_norm(k)=sum(sum((p-fid.p).^2,1),2);

    
	%%% Test minimum entropy

    e_k=y-Ga*xa_k_kmin;
    S_k=Ga*P_k_kmin*Ga.'+R;
    Term1=log(det(S_k))*length(y);
    Term2=0;

    for j=1:length(y)
        Term2=Term2+e_k(:,j).'/S_k*e_k(:,j);
    end

   objf_1(k)=Term1;
   objf_2(k)=Term2;
   
end

close all

figure; hold on; grid on;
plot(e_norm,p_norm,'-ob');
for k=1:length(Qp_val)
text(e_norm(k),p_norm(k),['Q_p=' num2str(Qp_val(k),'%0.0e')]);
end
xlabel('Error norm $||y-G \hat{x}-J \hat{p}||$','Interpreter','latex');
ylabel('Force norm $||\hat{p}||$','Interpreter','latex');

xlog; ylog;

plotopt=struct;
plotopt.LineStyleSet={'-' '--' '--' ':'}

ind=find(Qp_val==100);

[~,ind]=min(p_err_norm);

Qp_opt=Qp_val(ind)
% ind=ind-3

plotTime(fid.t,fid.p,p_all{ind},p_all{ind-2},p_all{ind+2},plotopt);

%% Test 

% close all

figure; hold on; grid on;
plot(Qp_val,p_err_norm,'-dk');
xlabel('Random walk $Q_p$','Interpreter','latex');
ylabel('Error norm $||p-\hat{p}||$','Interpreter','latex');

xlog; ylog;

figure; hold on; grid on;
plot(Qp_val,objf_1+objf_2,'-dk','DisplayName','Total sum');
plot(Qp_val,objf_1,'-ob','DisplayName','Model complexity');
plot(Qp_val,objf_2,'-xr','DisplayName','Data fit error');

xlabel('Random walk $Q_p$','Interpreter','latex');
ylabel('Negative LogLikelihood');
legend show

xlog

tilefigs



