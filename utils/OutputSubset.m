function [y_out] = OutputSubset(a_cell_all,a_cell_out,y_all);

%%

if size(y_all,1)~=length(a_cell_all);
    error('Length of source cell does not match signals in y');
end

if ischar(a_cell_out)
	a_cell_out={a_cell_out}
end

[dummy, ind_out1,ind_out2]=intersect(a_cell_all,a_cell_out,'stable');

[~,index_sorted]=sort(ind_out2);

ind_out=ind_out1(index_sorted); % the order same as in the requested a_cell_out
y_out=y_all(ind_out,:);