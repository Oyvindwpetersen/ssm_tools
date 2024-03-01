function [A B C D F]=ssmod_c2d(Ac,Bc,Cc,Dc,dt,varargin)

%% Continuous to discrete state-space model
%
% Inputs:
% Ac: state matrix
% Bc: input matrix
% Cc: output matrix
% Dc: direct transmission matrix
% dt: time discretization
%
% Outputs:
% A: state matrix
% B: input matrix
% C: output matrix
% D: direct transmission matrix
% F: matrix for first-order-hold
%
%% Input

p=inputParser;
addParameter(p,'fast',true,@islogical)
parse(p,varargin{:});
fast = p.Results.fast;

%% Discretize

A=[];
B=[];
F=[];

if ~isempty(Ac) & fast==true
    A=expmnorm(dt*Ac);
elseif ~isempty(Ac) & fast==false
    A=expm(dt*Ac);
end

if ~isempty(Bc)
    Ac_inv=eye(size(Ac)) / Ac;
    B=[A-eye(size(A))]*(Ac_inv*Bc);
    F=Ac_inv*(B-Bc*dt)*dt^-1;
end

C=[];
D=[];

if ~isempty(Cc)
    C=Cc;
end

if ~isempty(Dc)
    D=Dc;
end

