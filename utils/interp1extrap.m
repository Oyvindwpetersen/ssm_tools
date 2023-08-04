function Mi=interp1extrap(z,M,zi,method)
%%

if nargin==3
    method='linear';
end

if size(z,2)>size(z,1)
    z=z.';
end

if size(M,1)==size(z,2)
M=M.';
elseif size(M,2)==size(z,1)
M=M.';
end

if size(M,1)==size(M,2)
    warning('');
end

[z_sort,ind_sort]=sort(z);
M_sort=M(ind_sort,:);

%%

[min_z,indexMin]=min(z_sort);
[min_zi,~]=min(zi);

if min_zi<min_z
    z_expand_low=min_zi;
    M_expand_low=M_sort(indexMin,:);
else
    z_expand_low=[];
    M_expand_low=[];
end

[max_z,indexMax]=max(z_sort);
[max_zi,~]=max(zi);

if max_zi>max_z
    z_expand_up=max_zi;
    M_expand_up=M_sort(indexMax,:);
else
    z_expand_up=[];
    M_expand_up=[];
end

%%

M_expand=[ M_expand_low ; M_sort ; M_expand_up ];
z_expand=[ z_expand_low ; z_sort ; z_expand_up ];

Mi=interp1(z_expand,M_expand,zi,method,NaN);

% Mi=interp1(z,M,zi,method,NaN);

% indexReplace=find(isnan(Mi));
% 
% [min_z,indexMin]=min(z);
% [max_z,indexMax]=max(z);
% 
% indexBelow=indexReplace(zi(indexReplace)<min_z);
% indexAbove=indexReplace(zi(indexReplace)>max_z);
% 
% Mi(indexBelow)=M(indexMin);
% Mi(indexAbove)=M(indexMax);



end

