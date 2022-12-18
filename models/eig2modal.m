function [omega,xi]=eig2modal(lambda)

%% Convert eigenvalue to modal properties
%
% Inputs:
% lambda: eigenvalues, either array or cell with array
%
% Outputs:
% omega: undamped natural frequency
% xi: damping ratio
%
%% If A is a cell with many A-matrices, calculate for all


if ~iscell(lambda)
    omega=abs(lambda); xi=-real(lambda)./omega;
elseif iscell(lambda)
	for k=1:length(lambda)
		omega{k}=abs(lambda{k}); xi{k}=-real(lambda{k})./omega{k};
	end
end