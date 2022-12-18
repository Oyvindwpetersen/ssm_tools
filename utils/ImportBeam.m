function mod=ImportBeam(L,E,I,m,N_el,varargin)

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
addParameter(p,'Nm',[],@isnumeric) %springsupport
addParameter(p,'k0',50,@isnumeric)
addParameter(p,'k_end',500,@isnumeric)
parse(p,varargin{:})

type=p.Results.type;
Nm=p.Results.Nm;
k0=p.Results.k0;
k_end=p.Results.k_end;

%%

if strcmpi(type,'simplysupport')
    [Kg,Mg,f,phi,phi_full,K,M,K_red,M_red,S_red]=SimplySupportModal(L,E,I,m,N_el);
elseif strcmpi(type,'springsupport')
    [Kg,Mg,f,phi,phi_full,K,M,K_red,M_red,S_red]=SpringSupportModal(L,E,I,m,N_el,Nm,k0,k_end);
end

mod=struct();

mod.Kg=Kg;
mod.Mg=Mg;
mod.f=f;
mod.phi_red=phi;
mod.phi=phi_full;
mod.K=K;
mod.M=M;
mod.K_red=K_red;
mod.M_red=M_red;
mod.S_red=S_red;

mod.Omega=diag(mod.f)*2*pi;

dofs={'U' 'UR'};
number=1:(N_el+1);

mod.doflabel=DofLabel(dofs,number);

mod.nm=Nm;

%% Select limited modes

if ~isempty(Nm)
    
    range_keep=[1:Nm];
    
    mod.Kg=mod.Kg(range_keep,range_keep);
    mod.Mg=mod.Mg(range_keep,range_keep);
    mod.Omega=mod.Omega(range_keep,range_keep);
    mod.f=f(range_keep);
    mod.phi_red=phi(:,range_keep);
    mod.phi=phi_full(:,range_keep);

end