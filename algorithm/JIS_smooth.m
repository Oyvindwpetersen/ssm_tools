function [x_smooth p_smooth P_x_ss P_p_ss P_xp_ss P_px_ss]=JIS_smooth(A,B,G,J,y,R,Q,S,x0,P01,L,varargin)

%% Joint input and state estimation with smoothing
%
% Model
% x(k+1)=A*x(k)+B*p(k)+w(k);
% y=G*x(k)+J*p(k)+v(k);
%
% Inputs:
% A: state matrix
% B: input matrix
% G: output matrix
% J: direct transmission matrix
% y: output vector
% R: output noise covariance
% Q: state noise covariance
% S: mixed noise covariance
% x0: initial state estimate
% P01: initial state error covariance
% L: number of lags
%
% Outputs:
% x_smooth: smoothed state estimate
% p_smooth: smoothed input estimate
% P_x_ss: smoothed state error covariance
% P_p_ss: smoothed input error covariance
% P_xp_ss: mixed state/input error covariance
% P_px_ss: mixed state/input error covariance

%% Parse inputs

p=inputParser;
% p.KeepUnmatched=true;

addParameter(p,'text',true,@islogical)
addParameter(p,'plot',true,@islogical)
addParameter(p,'convtol',1e-6,@isnumeric)
addParameter(p,'minsteps',25,@isnumeric)
addParameter(p,'maxsteps',20e3,@isnumeric)

parse(p,varargin{1:end});

do_text=p.Results.text;
do_plot=p.Results.plot;
convtol=p.Results.convtol;
minsteps=p.Results.minsteps;
maxsteps=p.Results.maxsteps;

%%

ns=size(A,1);
% nm=size(A,1)/2;
nt=length(y);
nd=size(G,1);
np=size(B,2);
% convtol=1e-4;
minsteps=max(minsteps,L);


%% Error handling

if L==0
    error('L must be larger than zero (minimum one)');
end


%% Matrices

%% Extended observability, Toeplitz, N-matrix

tstart=tic;

Obs_L=extendedobs(A,G,L);

H_L=toeplitz_block(A,B,G,J,L);

N_L=N_matrix(A,G,L);

%% Union matrices

Id1_union=[eye(ns) zeros(ns,ns*(L-1))]; Id1_union=sparse(Id1_union);
Id2_union=[eye(np) zeros(np,np*L)]; Id2_union=sparse(Id2_union);
B_union=[B  zeros(ns,np*L)]; B_union=sparse(B_union);

%% Extended noise matrices

Q_L_i=cell(1,L);
R_Lplus_i=cell(1,L);
S_L_i=cell(1,L);
S_L_minusi=cell(1,L);

for i=1:L

    [Q_L_i{i},R_Lplus_i{i},S_L_i{i}]=extendednoisecov(Q,R,S,L,i);
    [~,~,S_L_minusi{i}]=extendednoisecov(Q,R,S,L,-i);

end

[Q_L_nil,R_Lplus_nil,S_L_nil]=extendednoisecov(Q,R,S,L,0);

Rbar_k_precalc=R_Lplus_nil+N_L*Q_L_nil*N_L.'+N_L*S_L_nil+S_L_nil.'*N_L.';
Sbar_k_precalc=Id1_union*Q_L_nil*N_L.'+Id1_union*S_L_nil;

telapsed=toc(tstart);
if do_text
    disp(['Base matrices calculated in ' sprintf('%2.1f', telapsed) ' s']);
end
% return

clear Pp_save M_L_k_save  Rbar_k_cond

%% Calculate steady state

tstart=tic;

P_x_k_kLminus=P01; %initial error covariance for state
K_L_k=cell(1,maxsteps);
CommonSumTerm_k=cell(1,maxsteps);

convreached=false; clear r_traceP;
Pp_prev=nan(np); Pp_current=nan(np);

% for k=0:(nt-1)
n=0; k=-1;
while convreached==false

    k=k+1; %time index, k=0,1,2,3...
    n=n+1; %matlab matrix index, n=1,2,3,4...

    %%% Input estimation

    %Pxw and Pxv
    Sum_loop_xw=zeros(ns,L*ns);
    Sum_loop_xv=zeros(ns,(L+1)*nd);
    Product_loop_i=eye(ns);

    if k>=1
        CommonSumTerm_k{n-1}=(Id1_union-K_L_k{n-1}*N_L);
    end

    for i=1:min(k,L)

        %     Product_loop_i=eye(ns);
        %     for j=1:(i-1)
        %         Product_loop_i=Product_loop_i*(A-K_L_k{n-j}*Obs_L);
        %     end

        if i>=2
            Product_loop_i=Product_loop_i*(A-K_L_k{n-(i-1)}*Obs_L); %Recursive PI-product
        end

        %     CommonSumTerm=(Id1_union-K_L_k{n-i}*N_L);
        CommonSumTerm=CommonSumTerm_k{n-i};

        SumTerm_xw_i=CommonSumTerm*Q_L_i{i}-K_L_k{n-i}*S_L_minusi{i}.';
        SumTerm_xv_i=CommonSumTerm*S_L_i{i}-K_L_k{n-i}*R_Lplus_i{i};

        Sum_loop_xw=Sum_loop_xw+Product_loop_i*SumTerm_xw_i;
        Sum_loop_xv=Sum_loop_xv+Product_loop_i*SumTerm_xv_i;

    end

    P_xw_k=Sum_loop_xw;
    P_xv_k=Sum_loop_xv;

    TempTerm=N_L*P_xw_k.'*Obs_L.'+P_xv_k.'*Obs_L.';
    Rbar_k=Obs_L*P_x_k_kLminus*Obs_L.'+Rbar_k_precalc+TempTerm+TempTerm.';
    Rbar_k=forcesym(Rbar_k);

    % if k>1
    Rbar_k_inv=Rbar_k\eye(size(Rbar_k));
    % else
    % Rbar_k_inv=pinv(Rbar_k);
    % end
    Rbar_k_inv=forcesym(Rbar_k_inv);

    HRbarH=H_L.'*Rbar_k_inv*H_L; HRbarH=forcesym(HRbarH);
    % HRbarH_pinv=pinv(HRbarH);
    HRbarH_pinv=HRbarH\eye(size(HRbarH));
    HRbarH_pinv=forcesym(HRbarH_pinv);

    M_L_k=Id2_union*HRbarH_pinv*H_L.'*Rbar_k_inv;

    P_p_k_kL=M_L_k*Rbar_k*M_L_k.'; P_p_k_kL=forcesym(P_p_k_kL);

    %%% State estimation

    Sbar_k=A*P_x_k_kLminus*Obs_L.'+Sbar_k_precalc+A*P_xw_k*N_L.'+Id1_union*P_xw_k.'*Obs_L.'+A*P_xv_k;

    TempTerm=Id1_union*P_xw_k.'*A.';
    Tbar_k=A*P_x_k_kLminus*A.'+Q+TempTerm+TempTerm.';
    Tbar_k=forcesym(Tbar_k);

    K_L_k{n}=Sbar_k*Rbar_k_inv - (Sbar_k*Rbar_k_inv*H_L-B_union)*HRbarH_pinv*H_L.'*Rbar_k_inv;

    P_x_kplus_kL=K_L_k{n}*Rbar_k*K_L_k{n}.'-K_L_k{n}*Sbar_k.'-Sbar_k*K_L_k{n}.'+Tbar_k; P_x_kplus_kL=forcesym(P_x_kplus_kL);

    %%% Next iteraton initialization

    P_x_k_kLminus=P_x_kplus_kL;

    %%% MISC SAVE

    % Pp_save(1,n)=P_p_k_kL;
    % M_L_k_save(:,n)=M_L_k;
    % Rbar_k_cond(n)=cond(Rbar_k);


    %%% Convergence checks

    Pp_prev=Pp_current;
    Pp_current=P_p_k_kL;

    r_traceP(n)=trace(abs(Pp_current-Pp_prev))./trace(Pp_current);

    if do_text & mod(k,10)==0  & k>0
        disp(['Ratio trace Pp ' num2str(r_traceP(n),'%0.3e')]);
    end

    if k>=minsteps & abs(r_traceP(n)) < convtol
        convreached=true;
        disp(['Trace convergence reached, k=' num2str(k)]);
        K_L_ss=K_L_k{n};
        M_L_ss=M_L_k;
        P_x_ss=P_x_kplus_kL;
        P_p_ss=P_p_k_kL;
        P_xp_ss=-(P_x_kplus_kL*Obs_L.'+P_xw_k*N_L.'+P_xv_k)*M_L_ss.';
        P_px_ss=P_xp_ss.';
    elseif k>= maxsteps
        disp(['Trace divergence reached, k=' num2str(k)]); break;
    end


end

telapsed=toc(tstart);

if do_text
    disp(['Steady state convergence calculated in ' sprintf('%2.1f', telapsed) ' s']);
end
%%

if strcmpi(do_plot,'yes')
    % close all;
    figure(); hold on; grid on;
    plot(r_traceP); set(gca,'YScale','log');
    xlabel('Steps');
    ylabel('Ratio trace Pp');
end

%% Calculate estimates

x_kplus_kL=zeros(ns,nt); x_kplus_kL(:,n)=x0;
phat_k_kL=zeros(np,nt);

n=0;
for k=0:(nt-1-L)
    n=n+1;

    % d_L_k=[d_L_k((nd+1):end) ; y(:,n+L) ];
    d_L_k=reshape(y(:,n:(n+L)),nd*(L+1),1);

    phat_k_kL(:,n)=M_L_ss*(d_L_k - Obs_L * x_kplus_kL(:,n));
    x_kplus_kL(:,n+1)=A*x_kplus_kL(:,n)+K_L_ss*(d_L_k-Obs_L*x_kplus_kL(:,n));

end

x_smooth=x_kplus_kL;
p_smooth=phat_k_kL;

end

%%

function [X_L_i]=blockdiag_L(X,L,i)

[n1,n2]=size(X);

% Number of matrices
n_matrices=L-abs(i);

% Fill the shifted diagonal (i'th below main diagonal) with ones
block_index=diag(ones(1,n_matrices),-i);
block_index(block_index==0)=2;

% Find the index of positions to fill
[i1,i2]=find(block_index==1);

% Fill the big matrix on the shifted diagonal
X_L_i=sparse(n1*L,n2*L);
for k=1:n_matrices
    range1=(i1(k)-1)*n1+[1:n1];
    range2=(i2(k)-1)*n2+[1:n2];

    X_L_i(range1,range2)=X;
end

end

%%

function checkdim(A,rows,cols)

[n1,n2]=size(A);

if n1~=rows & n2==cols
    error(['Row dimension of matrix wrong']);
elseif n1==rows & n2~=cols
    error(['Column dimension of matrix wrong']);
elseif n1~=rows & n2~=cols
    error(['Row and column dimension of matrix wrong']);
end

end

%%

function [Q_L,R_Lplus,S_L]=extendednoisecov(Q,R,S,L,i)


ns=size(Q,1);
nd=size(R,1);

Q_L=blockdiag_L(Q,L,i);

R_Lplus=blockdiag_L(R,L+1,i);

if i>=0
    S_temp=blockdiag_L(S,L,i);
    S_L=[S_temp zeros(L*ns,nd)];
elseif i<0
    S_temp=blockdiag_L(S,L,i+1);
    S_L=[zeros(L*ns,nd) S_temp];
end

if ( nnz(Q_L)/numel(Q_L) )<1e-4
    Q_L=sparse(Q_L);
end

if ( nnz(R_Lplus)/numel(R_Lplus) )<1e-4
    R_Lplus=sparse(R_Lplus);
end

if ( nnz(S_L)/numel(S_L) )<1e-4
    S_L=sparse(S_L);
end

checkdim(Q_L,L*ns,L*ns);

checkdim(R_Lplus,(L+1)*nd,(L+1)*nd);

checkdim(S_L,L*ns,(L+1)*nd);

end

%%
function Obs_L=extendedobs(A,G,L)

ns=size(A,1);
nd=size(G,1);

Obs_L=zeros((L+1)*nd,ns);

index_range=[1:nd];
for r=1:(L+1)
    Obs_L(index_range,:)=G*A^(r-1);
    index_range=index_range+nd;
end


end

%%

function N_L=N_matrix(A,G,L)

ns=size(A,1);
nd=size(G,1);

N0=zeros(nd,ns);
N1=[zeros(nd,ns) ; G];

N_k=cell(1,L);
N_k{1}=N1;

for k=2:L

    Obs_kminus=extendedobs(A,G,k-1);
    N_k{k}=zeros((k+1)*nd,k*ns);

    [i1,i2]=size(N_k{k-1});
    [j1,j2]=size(Obs_kminus);

    N_k{k}((end-i1+1):end,:)=[Obs_kminus N_k{k-1}];
end

if L==0
    N_L=N0;
elseif L==1
    N_L=N1;
elseif L>=2
    N_L=N_k{L};
end

end

%%

function H=toeplitz_block(A,B,G,J,L)

C=cell(1,L);

C{1}=J;
for r=2:(L+1)
    C{r}=G*A^(r-2)*B;
end

n=length(C);

[s1,s2]=size(J);
C{n+1}=zeros(s1,s2);

toeplitz_index=toeplitz(1:n);

for i=1:(n-1)
    for j=(i+1):n
        toeplitz_index(i,j)=n+1;
    end
end

C_stack=C(toeplitz_index);
H = cell2mat(C_stack);

end






