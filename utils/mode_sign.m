function [Phi,signf]=mode_sign(Phi,Phi_ref)

%% Change sign of mode shapes to match reference modes
%
% Inputs:
% Phi: matrix with modes shapes as columns
% Phi_ref: matrix with modes shapes as columns
%
% Outputs:
% Phi: matrix with modes shapes as columns (corrected)
% signs: vector with sign changes
%

%%

for k=1:size(Phi,2)

    c=corrcoef(Phi(:,k),Phi_ref(:,k));

    signs(k)=1;
    if c(2,1)<0
        Phi(:,k)=-Phi(:,k);
        signs(k)=-1;
    end

end
