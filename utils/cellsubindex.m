function index=cellsubindex(query_cell,parent_cell,old_or_new)

%% Function to find index of elements query_cell relative to elements in parent_cell
% Useful for finding selected DOFs among a list of all DOFs

% Example:
% parent_cell={'1_U1' '2_U1' '3_U1' '4_U1' '5_U1'}
% query_cell={'2_U1' '5_U1'}
% this will return
% index=[2 5]

% Inputs:
% query_cell: cell with strings of desired DOFs
% parent_cell: cell with strings of all DOFs
% old_or_new: type of algorithm, do not specify

% Outputs:
% index: vector with indices of query_cell in parent_cell

%% Check inputs

% Find algorithm to use
if nargin<3
    old_or_new='new';
    
    if length(query_cell)<20
        old_or_new='old';
    end
    
end

if ~iscell(query_cell)
    query_cell={query_cell};
end

if ~iscell(parent_cell)
    parent_cell={parent_cell};
end

index=zeros(1,length(query_cell));

%% Newer way, should be very fast for big cells

if strcmpi(old_or_new,'new')
    
    [truefalse, index] = ismember(query_cell,parent_cell);
    
    if ~any(truefalse) % no match, return
        index=nan(1,length(query_cell));
        warning(['***** query_cell not matched: ' 'none found']);
        return
    end
    
    index_missing=find(truefalse==0);
    
    if ~isempty(index_missing)
        index(index_missing)=NaN;
        for k=1:length(index_missing)
            warning(['***** query_cell not matched: ' query_cell{index_missing(k)}]);
        end
    end
    
end

%% Old way, search directly with strcmp

if strcmpi(old_or_new,'old')
    
    for k=1:length(query_cell)
        
        index_k=find(strcmp(parent_cell,query_cell{k}));
        
        if length(index_k)>1
            error(['***** Multiple matches in parent_cell for: ' query_cell{k}]);
        end
        
        if isempty(index_k)
            index(k)=NaN;
            warning(['***** query_cell not matched: ' query_cell{k}]);
        else
            index(k)=index_k;
        end
    end
    
end

%% Output

if size(index,1)>size(index,2)
    index=index.';
end

