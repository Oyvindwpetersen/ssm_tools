function [y_out,Gc,Jc]=PredictResponseModalForce(x,p,Mg,Cg,Kg,phi,dof_cell,doflabel,varargin)

%% Predict response (displacement or acceleration) from states and forces

%%
pi=inputParser;

addParameter(pi,'output','disp',@ischar)

parse(pi,varargin{:});

output = pi.Results.output;

%%

if isempty(dof_cell)
    error('No DOF');
end

if ~iscell(dof_cell)
    dof_cell={dof_cell};
end

if strcmpi(output,'disp')
    d_cell=dof_cell; a_cell={};
    
    [Sd,Sa,~]= DofSelection(d_cell,a_cell,{},doflabel);
    Gc=[Sd*phi , zeros(size(Sd*phi)) ];
    Jc=zeros(size(Gc,1),size(p,1));

    y_out=Gc*x;
    
elseif strcmpi(output,'acc')
    d_cell={}; a_cell=dof_cell;
    
    M=Mg;
    C=Cg;
    K=Kg;

    MinvK=M\K;
    MinvC=M\C;
    Minv=eye(size(M))/M;

    [Sd,Sa,~]= DofSelection(d_cell,a_cell,{},doflabel);
    Gc=[(Sd*phi-Sa*phi*MinvK) , -Sa*phi*MinvC];
    Jc=[Sa*phi * Minv];
    y_out=Gc*x+Jc*p;
end

%%







% if isempty(p); p=zeros(size(J,2),length(x)); end


% for k=1:ndof_est
%     
%    if strcmpi(output,'disp')
%         d_cell=dof_cell{k}; a_cell={};
% 	elseif strcmpi(output,'acc')
%         d_cell={}; a_cell=dof_cell{k};
%    end
% 
% 	[Sd,Sa,~]= generalSelection(d_cell,a_cell,{},doflabel);
%     
%     G=[Sd*phi-Sa*phi*Kg , -Sa*phi*Cg ];
%     J=[Sa*phi]; %*phi.'*Sp
%     
%     if isempty(p); p=zeros(size(J,2),length(x)); end
%     
% 	y_out(k,:)=G*x+J*p;
%     
% end