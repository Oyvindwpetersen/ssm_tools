function [ Mi ] = interp1z(z,M,zi,varargin)
%%interp1z Interpolates input data in similar fashion to interp1, but along z-axis (3. dimension)
%
% INPUT:    z:                  original z-axis
%           M:               	original matrix
%           zi:                 interpolated z-axis              
%                            
% OUTPUT:   Mi:                	interpolated matrix
%   
%
% Knut Andreas Kvaale (c) 2013
%

z=z(:);
zi=zi(:);

if nargin==3
    method='linear';
    extrapolate='extrap';
elseif nargin==5
    method=varargin{1};
    extrapolate=varargin{2};
end

[Lx,Ly,Lz] = size(M);
Lzi=length(zi);

Mmod=reshape(M,[],Lz).'; 

if strcmpi(extrapolate,'extrap_constant')
Mmodi = interp1extrap(z,Mmod,zi,method);
else
Mmodi = interp1(z,Mmod,zi,method,extrapolate);
end

Mi=reshape(Mmodi.',Lx,Ly,Lzi);


end

