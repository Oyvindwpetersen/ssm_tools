function H_RF=rf2tf(w_axis,a,d)

%% Transfer function for rational function

% Inputs:
% w_axis: frequency axis in [rad/s]
% a: scale factor vector
% d: poles (real)

% Outputs:
% H: transfer function

a_vec(:,1)=a;
d_vec(:,1)=d;

if size(w_axis,1)>size(w_axis,2)
    w_axis=w_axis.';
end

H_RF=a_vec./(1i*w_axis+d_vec);

H_RF=sum(H_RF,1);

H_RF=permute(H_RF,[1 3 2]);

