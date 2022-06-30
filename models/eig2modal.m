function [w,xi]=eig2modal(lambda)

%% Convert eigenvalue to modal properties

% Inputs:
% lambda: eigenvalues, either array or cell with array

% Outputs:
% w: undamped natural frequency
% xi: damping ratio

%% If A is a cell with many A-matrices, calculate for all


if ~iscell(lambda)
    w=abs(lambda); xi=-real(lambda)./w;
elseif iscell(lambda)
	for k=1:length(lambda)
		w{k}=abs(lambda{k}); xi{k}=-real(lambda{k})./w{k};
	end
end