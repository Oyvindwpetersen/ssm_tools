function [Sd,Sa,Sp]=DofSelection(d_cell,a_cell,p_cell,doflabel)

%% Create selection matrices

%%

if ~iscell(d_cell)
    d_cell={d_cell};
end

if ~iscell(a_cell)
    a_cell={a_cell};
end

if ~iscell(p_cell)
    p_cell={p_cell};
end

if length(unique(d_cell))~=length(d_cell)
    error('non-unique entries in d_cell')
end

if length(unique(a_cell))~=length(a_cell)
    error('non-unique entries in a_cell')
end

if length(unique(p_cell))~=length(p_cell)
    error('non-unique entries in p_cell')
end

ndof=length(doflabel);

dm=EstablishSelection(d_cell,doflabel);
am=EstablishSelection(a_cell,doflabel);
pm=EstablishSelection(p_cell,doflabel);

n_acc=size(am,2);
n_disp=size(dm,2);
n_force=size(pm,2);

ind_a=[1:n_acc];
ind_d=n_acc+[1:n_disp];
ind_p=[1:n_force];

np=n_force;
nd=n_acc+n_disp;

Sd=sparse(nd,ndof);
Sa=sparse(nd,ndof);
Sp=sparse(ndof,np);

if n_disp>0
Sd(ind_d',1:ndof)=dm';
end   

if n_acc>0
Sa(ind_a',1:ndof)=am';
end   

if np>0
Sp(1:ndof,ind_p)=pm;
end   

end

%%

function S_matrix=EstablishSelection(dof_cell,doflabel)

    if ~iscell(dof_cell)
        dof_cell={dof_cell};
    end
    
    if isempty(dof_cell)
        S_matrix=[]; return
    else
        [dof_index] = cellsubindex(dof_cell,doflabel);
    end
    
    s=ones(size(dof_index));
    S_matrix=sparse(dof_index,1:length(dof_index),s,length(doflabel),length(dof_index));

	r=nnz(S_matrix)/numel(S_matrix);
	if r>0.01
		S_matrix=full(S_matrix);
	end

end