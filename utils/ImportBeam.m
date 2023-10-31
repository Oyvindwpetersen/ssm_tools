function mod=importbeam(L,E,I,m,N_el,varargin)

%%

% Inputs: 
% 
% 

% Outputs: 
% 

%% Input handling

p=inputParser;
p.KeepUnmatched=true;

addParameter(p,'type','simplysupport',@ischar) %springsupport
addParameter(p,'nm',[],@isnumeric) %springsupport
addParameter(p,'k_end',0,@isnumeric)
parse(p,varargin{:})

type=p.Results.type;
nm=p.Results.nm;
k_end=p.Results.k_end;

%%

if strcmpi(type,'simplysupport')
    bc='support';
elseif strcmpi(type,'springsupport')
    bc='spring';
end

mod=struct();

[mod.Kg,mod.Mg,mod.f,mod.phi,mod.K,mod.M,mod.phi_red,mod.K_red,mod.M_red,mod.S_red,mod.S_red_inv,mod.nodecoord,mod.doflabel]= ...
simplysupportbeam(L,E,I,m,N_el,'bc',bc,'k_end',k_end);

mod.doflabel_red=mod.doflabel;
mod.doflabel_red(end-1)=[];
mod.doflabel_red(1)=[];


mod.Omega=diag(mod.f)*2*pi;

mod.nm=nm;

%% Select limited modes

if ~isempty(nm)
    
    range_keep=[1:nm];
    
    mod.Kg=mod.Kg(range_keep,range_keep);
    mod.Mg=mod.Mg(range_keep,range_keep);
    mod.Omega=mod.Omega(range_keep,range_keep);
    mod.f=mod.f(range_keep);
    mod.phi_red=mod.phi_red(:,range_keep);
    mod.phi=mod.phi(:,range_keep);

end