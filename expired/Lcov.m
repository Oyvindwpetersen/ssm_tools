%%

L_list=[0 1 2 3 4 5];
tic
parfor k=1:length(L_list)
    
L=L_list(k);

if L==0
[x_filt p_filt P_x_ss P_p_ss] = JIS_trunc_ss(fid.A,fid.B,fid.G,fid.J,fid.y,fid.x0,fid.R,fid.Q,fid.S,fid.P01);
else
[x_smooth p_smooth P_x_ss P_p_ss P_xp_ss P_px_ss ]=JIS_smooth(fid.A,fid.B,fid.G,fid.J,fid.y,fid.R,fid.Q,fid.S,fid.x0,fid.P01,L);
end
Pp{k}=diag(P_p_ss);
Px{k}=diag(P_x_ss);
end
toc

%%
close all

figure(); hold on;
for k=1:length(L_list)
plot(Pp{k}.^0.5./Pp{1}.^0.5,'-x'); %
end

set(gca,'YScale','log');

% xticks([1:length(fid.p_eq)]);
% xticklabels(ssmod.p_eq); xtickangle(90)


figure(); hold on;
for k=1:length(L_list)
plot(Px{k}.^0.5./Px{1}.^0.5,'-x');
end

set(gca,'YScale','log');

% xticks([1:length(ssmod.p_eq)]);
% xticklabels(ssmod.p_eq); xtickangle(90)


%%


close all
ssmod.L=30
H_L=toeplitzBlockMatrix(ssmod.A,ssmod.B,ssmod.G,ssmod.J,ssmod.L);
H_Lminus=toeplitzBlockMatrix(ssmod.A,ssmod.B,ssmod.G,ssmod.J,ssmod.L-1);
rank(H_L)-rank(H_Lminus)

plotSVD(H_L)