function H_RF=rf2tf(omega_axis,a,d)

%% Transfer function for rational function
%
% Inputs:
% omega_axis: frequency axis in [rad/s]
% a: scale factor vector
% d: poles (real)
%
% Outputs:
% H: transfer function

%%

a_vec(:,1)=a;
d_vec(:,1)=d;

if size(omega_axis,1)>size(omega_axis,2)
    omega_axis=omega_axis.';
end

H_RF=a_vec./(1i*omega_axis+d_vec);

H_RF=sum(H_RF,1);

H_RF=permute(H_RF,[1 3 2]);

