function [O]=ssmod_obs(Ac,Gc)

%% State space observability matrix
%
% Inputs:
% Ac: state matrix (cont)
% Gc: output matrix (cont)
%
% Outputs:
% O: observability matrix
%
%% 

O=[];
for k=1:size(Ac,1)
   O=[O;Gc*Ac^(k-1)]; 
end