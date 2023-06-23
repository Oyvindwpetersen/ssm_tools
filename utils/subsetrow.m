function a=subsetrow(A,a_label,A_label)
%% Acquire row subset of matrix A by filtering from labels
%
% Inputs:
% A: (n*m) matrix or cell
% a_label: cell with n_subset strings or numeric vector with labels of desired rows
% A_label: cell with n strings or numeric vector with labels of all rows
%
% Outputs:
% a: (n_subset*m) matrix with n_subset selected rows from A
%
%% Check

if isempty(a_label); a=[]; return; end

if length(A_label)~=size(A,1)
    error('Length of A_label must equal the number of rows in A')
end

if ischar(a_label)
    a_label={a_label};
end

if ischar(A_label)
    A_label={A_label};
end

%% Numeric labels

if isnumeric(a_label)
    
    ind_row=[];
    
    for k=1:length(a_label)
        
        ind_match=find(a_label(k)==A_label);
        
        if isempty(ind_match)
            error(['No match for label ' num2str(a_label(k)) ]);
        end
        
        if length(ind_match)>1
            error(['Multiple matches '  '(' num2str(length(ind_match)) ')' ' for label ' num2str(a_label(k)) ]);
        end
        
        ind_row(k)=ind_match;
        
    end
    
    a=A(ind_row,:);
    
end

%% Cell string labels

if iscell(a_label)
    
    if isempty(A_label{1}); a=[]; return; end
    
    ind_row=cellsubindex(a_label,A_label);
    a=A(ind_row,:);
    
end

end

