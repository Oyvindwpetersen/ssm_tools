function [ind dist]=indexargmin(x,x_val)

%% Find index of x that most closely matches x_val

% Inputs:
% x: vector with all values
% x_val: vector with desired values

% Outputs:
% ind: closest index
% dist: distances/residuals


%%

x_val=double(x_val);

for k=1:length(x_val)
[dist(k),ind(k)]=min(abs(x-x_val(k)));
end

