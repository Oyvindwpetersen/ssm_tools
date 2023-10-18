function [Hpy Hx0y Hx1y]=JIS_tf(A,B,G,J,Q,R,S,dt,omega_axis,varargin)
%% Transfer function for steady state operation of joint input and state estimator
%
% Model
% x(k+1)=A*x(k)+B*p(k)+w(k)
% y(k)=G*x(k)+J*p(k)+v(k)
%
% Inputs:
% A: state matrix
% B: input matrix
% G: output matrix
% J: direct transmission matrix
% Q: state noise covariance
% R: output noise covariance
% S: mixed noise covariance
% dt: time step
% omega_axis: frequency axis
%
% Outputs:
% Hpy: matrix with TF, output-to-force estimate
% Hx0y: matrix with TF, output-to-filter estimate
% Hx1y: matrix with TF, output-to-prediction estimate
%

%% Parse inputs

p=inputParser;
addParameter(p,'showtext',true,@islogical)
addParameter(p,'dispconv',true,@islogical)
addParameter(p,'trunc',false,@islogical)
addParameter(p,'convtol',1e-6,@isnumeric)

parse(p,varargin{:});

showtext=p.Results.showtext;
dispconv=p.Results.dispconv;
trunc=p.Results.trunc;
convtol=p.Results.convtol;

%%

% warning('Not updated yet');

ns=size(A,1);
ny=size(G,1);
np=size(B,2);

y_dummy=nan(ny,1);
x0=[];
P0=[];

[~,~,~,~,M_ss,K_ss,Kbar_ss] = JIS_ss(A,B,G,J,y_dummy,x0,Q,R,S,P0,'showtext',showtext,'dispconv',dispconv,'trunc',trunc,'convtol',convtol);

%% Old

% 
% M2=[M_ss ; L_ss ; zeros(ns,ny) ];
% 
% for k=1:length(omega_axis)
% 
%     M1=[ eye(np) zeros(np,ns) M_ss*G ;
%         L_ss*J eye(ns) (-eye(ns)+L_ss*G);
%         B A -exp(1i.*omega_axis(k)*dt)*eye(ns) ];
% 
%     M3=M1\M2;
% 
%     M4(:,:,k)=M3; %M4 is [Hpd;Hx0d;Hx1d]
% 
% end

%%

M2=[M_ss ; K_ss ; Kbar_ss ];

for k=1:length(omega_axis)

    z=exp(1i.*omega_axis(k)*dt);

    M1=[ eye(np) zeros(np,ns) M_ss*G ;
        K_ss*J eye(ns) (-eye(ns)+K_ss*G);
        zeros(ns,np) zeros(ns,ns) z*eye(ns)-A+Kbar_ss*G ];

    M3=M1\M2;

    M4(:,:,k)=M3; %M4 is [Hpd;Hx0d;Hx1d]

end

Hpy=M4(1:np,:,:);
Hx0y=M4(np+[1:ns],:,:);
Hx1y=M4(np+ns+[1:ns],:,:);

