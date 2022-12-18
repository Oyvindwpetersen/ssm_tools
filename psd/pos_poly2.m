function [r,re,im]=pos_poly2(c_even)

%% Find real and imaginary part of roots of an even polynomial
%
% Inputs:
% c_even: polynomial coeffcients for even terms, c_even=[c(n) c(n-2) ... c(2) c(0)];
%
% Outputs:
% r: all roots
% re: real parts
% im: imaginary parts
%
%% Order of polynomial

% Add zeros for odd powers
c_all=reshape([c_even ; zeros(size(c_even))],1,[]); c_all(end)=[];

r=roots(c_all);

[~,ind_sort]=sort(real(r),'descend');

r_sort=r(ind_sort); r_sort=r_sort(1:2:length(r_sort)/2);

re=real(r_sort);
im=imag(r_sort);

order_c=length(c_all)-1;

if any(order_c==[2 6 10 14 18 22 26 30])
    
    re_take_a_closer_look=re(end);
    
    if abs(re_take_a_closer_look)>1e-12
        c_all
        r
        re
        re_1=re(1)
        re_end=re(end)
%         re_1=re(1)
        error('This number should be zero (purely imaginary pole). Please check.');
    end
    
    re(end)=[];
    
end

