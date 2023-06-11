function [r,c_even,c_all]=pos_poly(R,I,forcezero)

%% Create positive only polynomial from real and imaginary part of roots
% The polynomial will only have even terms, y=c0+c2*x^2+c4*x^4...
%
% Inputs:
% R: vector, real part of roots
% I: vector, imaginary part of roots
% forcezero: true/false, if true then two zero roots are added so that y(x=0)=0
%
% Outputs:
% r: all roots
% c_even: polynomial coeffcients for even terms, c_even=[c(k) c(k-2) ... c(2) c(0)];
% c_all: polynomial coeffcients for all terms, c=[c(k) 0 c(k-2) 0 ... 0 c(2) 0 c(0)];
%

%%

if nargin==2
    forcezero=false;
end

%% Order of polynomial

%Highest order of polynomial
order_poly=(length(R)+length(I))*2;

n_var=order_poly/2;

% Number of imaginary and real numbers (unique)
n_im=ceil(n_var/2);
n_re=n_var-n_im;

% Checks
if length(R)~=n_re
    n_re
    R
    error('Wrong size of real vector');
end

if length(I)~=n_im
    n_im
    I
    error('Wrong size of imaginary vector');
end
   
%%

% if order_poly==[2 6 10 14 ...], then there is two poles purely imaginary, the rest are in conjugate pairs of four (n_im=n_re+1)

% if order_poly==[4 8 12 16 ...], then the poles are in conjugate pairs of four (n_im=n_re)

% Examples:

% order=2
% a=[1 0 2 0 1.1]
% a(end)=0;
% r=roots(a)

% order=4
% a=[1 0 3 0 2 0 4];
% a(end)=0;
% r=roots(a)

% order=6
% a=[1 0 -3 0 5 0 -2 0 1];
% a(end)=0
% r=roots(a)

% order=8
% a=[1 0 -3 0 5 0 -2 0 3 0 1];
% a(end)=0
% r=roots(a)

%%

n_roots=order_poly;

r=zeros(n_roots,1);

% Produce conjugate pairs of four
for k=1:n_re
    
    range_add=(k-1)*4+[1:4];
    
    r_add=[...
        R(k)+1i*I(k) ; 
        R(k)-1i*I(k) ; 
        -R(k)+1i*I(k) ; 
        -R(k)-1i*I(k) ];
    
    r(range_add)=r_add;
    
end

% Add last root, imaginary only
if n_im>n_re
    
    range_add=[n_roots-1 n_roots];

    r_add=[...
        0+1i*I(end) ;
        0-1i*I(end) ];
    
    r(range_add)=r_add;
    
end

% Add two zero roots
if forcezero
    r=[r; 0 ; 0];
end

% Polynomial coefficients
c_all=poly(r);

% Odd coefficients should be zero
c_odd=c_all(2:2:end);
if any(abs(c_odd)./max(abs(c_all))>1e-12)
    
    c_odd
    c_all
    error('Coefficients for odd powers not zero');
    
end

% Check, coeff should be real coeff should be zero
if ~isreal(c_all) & any(imag(c_all)./max(abs(c_all))>1e-12)
	
    c_all
    error('Imaginary polynomial coeffcients, something wrong here');
    
end

% Set odd to zero
c_all(2:2:end)=0;
c_even=c_all(1:2:end);

