%%
close all
fid.L=2
H_L=toeplitzBlockMatrix(fid.A,fid.B,fid.G,fid.J,fid.L);
H_Lminus=toeplitzBlockMatrix(fid.A,fid.B,fid.G,fid.J,fid.L-1);
rank(H_L)-rank(H_Lminus)

% figure();

% fid.L=5


% 
plotSVD(H_L);
plotSVD(H_Lminus);

