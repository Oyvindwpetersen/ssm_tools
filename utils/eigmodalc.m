function [omega,xi,phi]=eigmodalc(K,C,M,normalize)

%% Find mode shapes
%
% Inputs:
% K: stiffness matrix
% C: damping matrix
% M: mass matrix
% normalize: 'mass'
%
% Outputs:
% omega: natural frequencies
% xi: damping ratios
% phi: mode shapes
%
%%

if nargin<4
    normalize='mass';
end

Ac=[zeros(size(K)) eye(size(K)) ; -M\K -M\C];

[lambda,psi,omega,xi]=eigA(Ac);

n=size(K,1);

xi=xi(1:2:end);
omega=omega(1:2:end);
phi=psi(1:n,1:2:end);

if strcmpi(normalize,'mass')

    [mpcw,rotrad,phir,mpc]=mpcweighted(phi);

    if all(mpcw>99.999)
        V=real(phir);
        for i=  1:n
            alfa=1/sqrt(V(:,i).'*M*V(:,i));
            phi(:,i)=V(:,i)*alfa;
        end

    else
        warning('Mode shapes not real, could not be mass normalized');
    end

    % Set max to positive value
    for k=1:size(phi,2)
        [~,ind_max]=max(abs(phi(:,k)));

        if phi(ind_max,k)<0
            phi(:,k)=-phi(:,k);
        end
    end

end

end

function [mpcw,rotrad,psir,mpc]=mpcweighted(psi)
%% Rotation of mode shapes for best alignment with real axis
%
% Inputs:
% psi: [ndof,nm] complex mode shapes
%
% Outputs:
% mpcw: [1,nm] weighted modal phase collinearity
% rotrad: [1,nm] rotation of mode shapes for best alignment real axis
% psir: [ndof,nm] complex rotated modes
% mpc: [1,nm] modal phase collinearity
%

%%

[ndof,nm]=size(psi);

% If empty, return
if isempty(psi)
    mpcw=NaN*ones(1,nm);
    mpc=NaN*ones(1,nm);
    rotrad=NaN*ones(1,nm);
    psir=NaN*ones(size(psi));
    return
end

% If real for all modes, then set to 100
if isreal(psi)
    mpcw=100*ones(1,nm);
    mpc=100*ones(1,nm);
    rotrad=0*ones(1,nm);
    psir=psi;
    return
end

% sxx=ones(1,ndof)*(real(psi).^2);
% syy=ones(1,ndof)*(imag(psi).^2);
% sxy=ones(1,ndof)*(real(psi).*imag(psi));


sxx=sum(real(psi).^2,1);
syy=sum(imag(psi).^2,1);
sxy=sum(real(psi).*imag(psi),1);

eta=(syy - sxx) ./ (2 * sxy);
% term0=sqrt((eta .* eta) + 1);
term0=sqrt((eta.^2) + 1);


beta=eta + (sign(sxy) .* term0);
term1=(sxx + syy) ./ 2;
term2=sxy .* term0;
ev1=term1 + term2;
ev2=term1 - term2;
eratio=(ev1 - ev2) ./ (ev1 + ev2);
tau=atan(beta);
mpcw0=eratio .* eratio;
mpcw=100*mpcw0;
rotrad0=tau + pi/2; % Rotation of mode shapes for best alignment with +/- 90 deg.
rotrad=rotrad0-pi/2; % Alignment with real axis

%% MPC

% (abs value of eta or not? Different practices exist in literature)

if nargout>3
    term_1=abs(eta)+sign(eta).*term0;
    term_2=atan(term_1);
    mpc=(sxx+1./eta.*sxy.*(2*(eta.*eta+1).*sin(term_2).*sin(term_2)-1))./(sxx+syy);
    mpc=mpc*100;
end

%% Catch single modes that might be NaN

if any(isnan(mpcw) | isinf(mpcw))

    idx_check=find(isnan(mpcw) | isinf(mpcw));
    for k=1:length(idx_check)

        % If real, then set to 100
        if isreal(psi(:,idx_check(k)))
            mpcw(idx_check(k))=100;
            mpc(idx_check(k))=100;
            rotrad(idx_check(k))=0;
            psir(:,idx_check(k))=psi(:,idx_check(k));
            return
        else
            error('Mode returned NaN, check here');
        end

    end
end

%% Alignment with real axis

% psir_old=conj(psi) .* repmat(exp(1i*rotrad),ndof,1);

if nargout>2
    psir=conj(psi) .* exp(1i*rotrad); % Without repmat (still works, matches number of columns)
end

end