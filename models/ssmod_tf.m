function H=ssmod_tf(A,B,C,D,w,dt,varargin)

%% Transfer function for state space model (input-output or input-state)
%
% Inputs:
% A: state matrix (cont or disc)
% B: input matrix (cont or disc)
% C: output matrix (cont or disc)
% D: direct transmission matrix (cont or disc)
% w: frequency vector in rad/s
% dt: time discretization, set to empty [] if cont
%
% Outputs:
% H: transfer function
%
%%

p=inputParser;
addParameter(p,'type','io',@ischar)
parse(p,varargin{:});
type=p.Results.type;

%%

if nargin==5 | isempty(dt)
    dt=[]; time_type='cont';
elseif nargin==6
	time_type='disc';
end

np=size(B,2);
ns=size(A,1);
ny=size(C,1);

%%

if strcmpi(type,'io')

        H=zeros(ny,np,length(w));

    if time_type=='cont'

        for j=1:length(w)
            H(:,:,j)=C / (-A+1i*w(j)*eye(ns)) * B+D;
        end

	elseif time_type=='disc'

        for j=1:length(w)
            H(:,:,j)=C / (-A+exp(1i*w(j)*dt)*eye(ns)) * B+D;
        end  

    end

end

%%

if strcmpi(type,'is')

    H=zeros(ns,np,length(w));

    if time_type=='cont'

        for j=1:length(w)
            H(:,:,j)=eye(ns) / (-A+1i*w(j)*eye(ns)) * B;
        end

    elseif time_type=='disc'

        for j=1:length(w)
            H(:,:,j)=eye(ns) / (-A+exp(1i*w(j)*dt)*eye(ns)) * B;
        end  

    end

end
