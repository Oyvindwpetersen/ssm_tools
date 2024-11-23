function [label_cell] = genlabels2(dofs,number,varargin)
%% Produce cell of labels of degrees of freedom
%
% Inputs: 
% dofs: dof labels, e.g. {'U1' 'U2' 'U3'}
% number: node numbers, e.g. [10:20]
%
% Outputs: 
% label_cell: labels, e.g.  {'10_U1' '10_U2' '10_U3' '11_U1' ... '20_U3' }  
%
%% Test data
%
%[dof_cell] = genlabels({'U1' 'U2' 'U3'},[ 410 412])
%
%% Input handling

% p=inputParser;
% p.KeepUnmatched=true;
% 
% addParameter(p,'postfix','',@ischar)
% addParameter(p,'prefix','',@ischar)
% parse(p,varargin{:})
% 
% postfix=p.Results.postfix;
% prefix=p.Results.prefix;

postfix='';
prefix='';

%%

% If empty, return empty
if isempty(dofs) | isempty(number)
    label_cell={''};
    return
end


%Put dof label in cell
if ischar(dofs)
    dofs={dofs};
end

%Special cases
if ~isempty(dofs) & length(dofs)>1
    % Ok
elseif isempty(dofs) 
        dofs={'U1' 'U2' 'U3' 'UR1' 'UR2' 'UR3'};
else
    if strcmpi(dofs(1),'all')
        dofs={'U1' 'U2' 'U3' 'UR1' 'UR2' 'UR3'};
    elseif strcmpi(dofs(1),'trans')
        dofs={'U1' 'U2' 'U3' };
    elseif strcmpi(dofs(1),'rot')
        dofs={'UR1' 'UR2' 'UR3'};
    end
end

%Create cell
n1=length(number);
n2=length(dofs);

label_cell=cell(n1*n2,1);
idx=[1:length(dofs)];
for k=1:n1
    num=sprintf(['%d'],number(k)); % from int2str

    fullLabels = cellfun(@(x)[[num '_'], x], dofs, 'UniformOutput', false);

    label_cell(idx)=fullLabels;

    idx=idx+length(dofs);
    % for j=1:n2
    %     i=i+1;
    %     label_cell{i}=[prefix num '_' dofs{j} postfix ];
    % end
end

