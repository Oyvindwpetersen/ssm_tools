function [p,q]=rf2poly(a,d)

%% Poly coefficients for rational function

% Inputs:
% a: scale factor vector
% d: poles (real)

% Outputs:
% p: polynomial coefficients
% q: polynomial coefficients

%%
% H(w)= sum_k a_k/(iw+d_k) = a_1/(iw+d_1) + a_2/(iw+d_2) + a_3/(iw+d_3)...
%
% Common denominator form:
% H=TOP/BOTTOM
% BOTTOM = (iw+d_1)(iw+d_2)(iw+d_3): product of all
% TOP = a_1*(iw+d_2)(iw+d_3)+a_2*(iw+d_1)(iw+d_3)+a_3*(iw+d_1)(iw+d_2): one removed
% H=polynomial_of_p / polynomial_of_q

%%
q=poly(-d);

% r=length(d)
% q=[q_r q_r-1 ... q0]

r=length(d);

clear p_component
for k=1:r

    d_component=d; d_component(k)=[];
    p_component{k}=poly(-d_component)*a(k);

end

p_component=cell2mat(p_component.');
p=sum(p_component,1);

%% test

% a=[5 20]
% d=[2 6];
% 
% w=[0:0.1:100];
% 
% H_poly=poly2tf(w,p,q);
% 
% H_RF=rf2tf(w,a,d);
% 
% 
% close all
% plotTransferFunction(w,H_poly,w,H_RF);