%% Estimate

close all

[x_jis p_jis Px_jis_ss Pp_jis_ss ] = JIS_trunc_ss(mod.A,mod.B,mod.G,mod.J,fid.y,fid.x0,fid.R,fid.Q,fid.S,fid.P01,'trunc','yes');

fid.L=60;
[x_smooth p_smooth Px_smooth_ss Pp_smooth_ss ]=JIS_smooth(mod.A,mod.B,mod.G,mod.J,fid.y,fid.R,fid.Q,fid.S,fid.x0,fid.P01,fid.L);


plotTime(fid.t,fid.p,p_jis,p_smooth,'legend',{'True' 'Filt' 'Smooth'});

plotFreq(fid.t,fid.p,p_jis,p_smooth,'xlim',[0 50],'legend',{'True' 'Filt' 'Smooth'});

plotTime(fid.t,fid.x(1:6,:),x_jis(1:6,:),x_smooth(1:6,:));

plotFreq(fid.t,fid.x(1:6,:),x_jis(1:6,:),x_smooth(1:6,:),'xlim',[0 50]);

tilefigs



%% Check statistics (best for long time series, band limited white noise input)

err_x_jis=std(fid.x-x_jis,0,2) %./ std(fid.x,0,2);
err_x_smooth=std(fid.x-x_smooth,0,2) %./ std(fid.x,0,2);

err_p_jis=std(fid.p-p_jis,0,2) %./ std(fid.p,0,2);
err_p_smooth=std(fid.p-p_smooth,0,2) %./ std(fid.p,0,2);

close all

% figure(); hold on; grid on; ylog;
% plot(err_x_jis,'ob','DisplayName','Error JIS');
% plot(err_x_smooth,'xr','DisplayName','Error Smooth');
% plot(diag(Px_jis_ss).^0.5,'*k','DisplayName','Error JIS (from filter)');
% plot(diag(Px_smooth_ss).^0.5,'dk','DisplayName','Error Smooth (from filter)');
% legend show
% 
% figure(); hold on; grid on; ylog;
% plot(err_p_jis,'ob','DisplayName','Error JIS');
% plot(err_p_smooth,'xr','DisplayName','Error Smooth');
% plot(diag(Pp_jis_ss).^0.5,'*k','DisplayName','Error JIS (from filter)');
% plot(diag(Pp_smooth_ss).^0.5,'dk','DisplayName','Error Smooth (from filter)');
% legend show


figure(); hold on; grid on; ylog;
title('Ratio variance state'); legend show;
plot(err_x_jis./diag(Px_jis_ss).^0.5,'ob','DisplayName',' JIS');
plot(err_x_smooth./diag(Px_smooth_ss).^0.5,'*r','DisplayName',' Smooth');


figure(); hold on; grid on; ylog;
title('Ratio variance'); legend show;
plot(err_p_jis./diag(Pp_jis_ss).^0.5,'ob','DisplayName',' JIS');
plot(err_p_smooth./diag(Pp_smooth_ss).^0.5,'*r','DisplayName',' Smooth');

% 
% figure(); hold on; grid on; ylog;
% plot(err_p_jis,'ob','DisplayName','Error JIS');
% plot(err_p_smooth,'xr','DisplayName','Error Smooth');
% plot(diag(Pp_jis_ss).^0.5,'*k','DisplayName','Error JIS (from filter)');
% plot(diag(Pp_smooth_ss).^0.5,'dk','DisplayName','Error Smooth (from filter)');
% legend show
% 
% 

tilefigs


%% Test lag

% L_mat=[1 3 5 10 20 30 40 ];
L_mat=[10:10:80];

for k=1:length(L_mat)
    
    [~,~,P_x_ss,P_p_ss ]=JIS_smooth(mod.A,mod.B,mod.G,mod.J,fid.y,fid.R,fid.Q,fid.S,fid.x0,fid.P01,L_mat(k));
    Px_all(k,:)=diag(P_x_ss);
    Pp_all(k,:)=diag(P_p_ss);
end

close all
figure(); hold on; grid on;
plot(L_mat,Px_all); ylog;
    
figure(); hold on; grid on;
plot(L_mat,Pp_all); ylog;
