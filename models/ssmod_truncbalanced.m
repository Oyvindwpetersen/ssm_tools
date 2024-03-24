function [Ac_t Bc_t Cc_t Dc_t T]=ssmod_truncbalanced(Ac,Bc,Cc,Dc,n,varargin)

%% Balanced trunction for state space model (input-output or input-state)
%
% See Unneland, Chapter 5.2
%
% Inputs:
% Ac: state matrix (cont)
% Bc: input matrix (cont)
% Cc: output matrix (cont)
% Dc: direct transmission matrix (cont)
% n: truncation order (dimension of truncated problem)
%
% Outputs:
% Ac_t: state matrix (cont) (truncated)
% Bc_t: input matrix (cont) (truncated)
% Cc_t: output matrix (cont) (truncated)
% Dc_t: direct transmission matrix (cont) (truncated)
%
%%

p=inputParser;
p.KeepUnmatched=true;
addParameter(p,'plot',true,@islogical)

parse(p,varargin{:})
plot_sv=p.Results.plot;

%%

ss_full=ss(Ac,Bc,Cc,Dc);

contGram = gram(ss_full,'c'); contGram=forcesym(contGram);
obsGram = gram(ss_full,'o'); obsGram=forcesym(obsGram);

M=contGram; N=obsGram; 
% N=N+diag(diag(N))*1e-9;
Lm = chol(M,'lower'); M_rec=Lm*Lm.';
Ln = chol(N,'lower'); N_rec=Ln*Ln.'; 

[U,S,V]=svd(Lm.'*Ln); SD=diag(S);

T=S.^0.5*U.' / Lm;
Tinv=eye(size(Ln)) / (Ln.') *V*S.^0.5;

Abar=T*Ac*Tinv;
Bbar=T*Bc;
Cbar=Cc*Tinv;
Dbar=Dc;

Ac_t=Abar(1:n,1:n);
Bc_t=Bbar(1:n,:);
Cc_t=Cbar(:,1:n);
Dc_t=Dbar;

if strcmpi(plot_sv,true)

SD_sum=sum(SD);

figure(); hold on; grid on;
plot(SD./SD_sum,'x');




SD_cumsum=cumsum(SD);
plot(SD_cumsum./SD_sum,'o');

ind=find(SD_cumsum./SD_sum>0.9,1); plot(ind,SD_cumsum(ind)./SD_sum,'dk');
ind=find(SD_cumsum./SD_sum>0.99,1); plot(ind,SD_cumsum(ind)./SD_sum,'dk');
ind=find(SD_cumsum./SD_sum>0.999,1); plot(ind,SD_cumsum(ind)./SD_sum,'dk');

title('Singulvar values (normalized)');
legend({'SV' 'Cummulative SV' '90%' '99%' '99.9%'});

set(gca,'Yscale','log');


figure(); hold on; grid on;
plot(SD./SD(1),'x'); ylog;

end