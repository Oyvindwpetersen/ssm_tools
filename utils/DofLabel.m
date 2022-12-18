function [label_cell] = DofLabel(dofs,number,varargin)
%%

% Inputs: 
% dofs: dof labels, e.g. {'U1' 'U2' 'U3'}
% number: node numbers, e.g. [10:20]

% Outputs: 
% label_cell: labels, e.g.  {'10_U1' '10_U2' ... '20_U3' }  

%% Test data

%[dof_cell] = getLabel({'U1' 'U2' 'U3'},[ 410 412])

%% Input handling

p=inputParser;
p.KeepUnmatched=true;

addParameter(p,'postfix','',@ischar)
addParameter(p,'prefix','',@ischar)
parse(p,varargin{:})

postfix=p.Results.postfix;
prefix=p.Results.prefix;

%%

if isempty(dofs) | isempty(number)
    label_cell={''};
    return
end

%Put dof number in cell
if isnumeric(number)
    number_cell=cell(1,length(number));
    for k=1:length(number)
        number_cell(k)={num2str(number(k))};
    end
elseif iscell(number)
    number_cell=number;
end

%Put dof label in cell
if ischar(dofs)
    dofs={dofs};
end

%Special cases
if isempty(dofs) 
    dofs={'U1' 'U2' 'U3' 'UR1' 'UR2' 'UR3'};
elseif strcmpi(dofs(1),'all')
    dofs={'U1' 'U2' 'U3' 'UR1' 'UR2' 'UR3'};
elseif strcmpi(dofs(1),'trans')
    dofs={'U1' 'U2' 'U3' };
elseif strcmpi(dofs(1),'rot')
    dofs={'UR1' 'UR2' 'UR3'};
end

%Create cell
n1=length(number_cell);
n2=length(dofs);


label_cell=cell(1,n1*n2);
i=0;
for k=1:n1
    for j=1:n2
        i=i+1;
        label_cell{1,i}=char([prefix char(number_cell(k)) '_' upper(char(dofs(j))) postfix ]);
    end
end
