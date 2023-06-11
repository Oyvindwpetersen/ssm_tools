function [Px Pp Pd] = variableSmoothingCov(A,B,G,J,y,R,Q,S,x0,P01,phi,omega,gamma,Sp,d_ref,a_ref,doflabel,L_list,varargin);

%%

% pi=inputParser;
% addParameter(pi,'force','disc',@ischar)

% parse(pi,varargin{:});
% forcetype = pi.Results.force;

% forcetype_check=[ strcmpi(forcetype,'disc') strcmpi(forcetype,'discrete') strcmpi(forcetype,'modal')];
% if sum(forcetype_check)==0
    % error();
% end


%
% [Sde,Sae,~]= generalSelection(d_ref,a_ref,{},doflabel);

% if strcmpi(forcetype,'disc') | strcmpi(forcetype,'discrete');
    % [~,~,Ge,Je]=statespaceModel(phi,omega,gamma,Sae,Sde,Sp,0);
% elseif strcmpi(forcetype,'modal')
    % [~,~,Ge,Je]=statespaceModelModalForce(phi,omega,gamma,Sae,Sde,Sp,0);
% end

% for k=1:length(L_list)
    
% [~,~,P_x_ss{k} P_p_ss{k} P_xp_ss{k} P_px_ss{k}]=JIS_smooth(A,B,G,J,y,R,Q,S,x0,P01,L_list(k));

% P_de{k}=[Ge Je]*[P_x_ss{k} P_xp_ss{k} ; P_px_ss{k} P_p_ss{k}]*[Ge.' ; Je.'];

% Px(:,k)=diag(P_x_ss{k});
% Pp(:,k)=diag(P_p_ss{k});
% Pd(:,k)=diag(P_de{k});

% end


% return
%

% close all;
% figure(); hold on; grid on;
% plot(Pp'); set(gca,'YScale','log');

% figure(); hold on; grid on;
% plot(Px'); set(gca,'YScale','log');

% figure(); hold on; grid on;
% plot(Pd'); set(gca,'YScale','log');

