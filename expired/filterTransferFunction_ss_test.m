
%%

ssmod=fid;
fid.P0=fid.P01
ssmod.ns=fid.nm*2;
ssmod.ny=fid.nd;
ssmod.dt=dt;

[fid.xfilt fid.pfilt,~,~,M_ss,L_ss] = JIS_trunc_ss(ssmod.A,ssmod.B,ssmod.G,ssmod.J,fid.y,fid.x0,fid.R,fid.Q,fid.S,fid.P0,'trunc','no');

% [fid.xfilt fid.pfilt,~,~,M_ss,L_ss] = JIS_trunc_ss(fid.A,fid.B,fid.jis.G,fid.jis.J,fid.jis.y,fid.x0,fid.jis.R,fid.Q,fid.jis.S,fid.P0,'trunc','no');
%%
clear Hpd M1 M2 M3 M4

M2=[M_ss ; L_ss ; zeros(ssmod.ns,ssmod.ny) ];
omega=[0:0.5:300];

for k=1:length(omega)
    
M1=[  eye(ssmod.np) zeros(ssmod.np,ssmod.ns) M_ss*ssmod.G ;
L_ss*ssmod.J eye(ssmod.ns) (-eye(ssmod.ns)+L_ss*ssmod.G);
ssmod.B ssmod.A -exp(1i.*omega(k)*ssmod.dt)*eye(ssmod.ns) ];

M3=eye(size(M1)) / M1 * M2; 

M4(:,:,k)=M3; %M4 is [Hpd;Hx0d;Hx1d]^T

Hpd(:,:,k)=M4(1:ssmod.nm,:,k); %M4 is [Hpd;Hx0d;Hx1d]^T

end

%%
% close all
figure(); makebig();
ha = tight_subplot( 1,5,[.05 .05],[.1 .1],[.1 .05]);

modeno=3

for jj=1:length(ha)
    
    axes(ha(jj)); hold on; grid on;
    set(gca,'YScale','log');

    plot(omega,abs(squeeze(Hpd(modeno,jj+0,:))));
    plot(omega,abs(squeeze(Hpd(modeno,jj+5,:))));

    axistight(gca,[0 0.05],'x','ylog');
%     ylabel(ssmod.a_o{jj});
    if jj==1; legend({'Acc' 'Disp'}); end
%     ylim([min(min(abs(Hpd(modeno,1:70,:)),[],3)) max(max(abs(Hpd(modeno,1:70,:)),[],3))]);
% 	ylim([10 4e3]);
end
%%

% %%
% [Sdd,omega_d]=estimateSpectrumWelch(fid.y,ssmod.Fs,'plot','no','unit','rad'); 
% 
% Hpd_int=interp1zbig(omega,Hpd,omega_d,'linear',0);
% Spp=mtimes3(Hpd_int,Sdd,Hpd_int,'nnh','mmx');
% 
% close all
% plotSpectrum(wave.w,wave.Spp_modal,omega_d,Spp,'xlim',[0 8],'log','no','component',[1:8]);
% plotSpectrum(wave.w,wave.Spp_modal,omega_d,Spp,'xlim',[0 8],'log','no','component',[9:15]);
% 
