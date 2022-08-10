function X_newton=dare_sparse_newton(A,H,Q,R,S,varargin)

%%

% Algorithm 1 from 
% On the numerical solution of large-scale sparse discrete-time Riccati equations
% Peter Benner & Heike Faßbender 
% Advances in Computational Mathematics volume 35, Article number: 119 (2011) 

% Solves X in the equation
% R(X)=0=Q+AXA'-X-(AXH'+S')*inv(R+H*X*H.')*(AXH'+S');
% Note that S might be transposed compared to other definitions (re-transpose automatic below) 

% A and H are sparse and large matrices

%% Parse inputs

p=inputParser;

addParameter(p,'max_iter',20,@isnumeric)
addParameter(p,'tol',1e-10,@isnumeric)
addParameter(p,'X0',[],@isnumeric)
addParameter(p,'compare',false)

parse(p,varargin{1:end});

max_iter=p.Results.max_iter;
tol=p.Results.tol;
X0=p.Results.X0;
DoCompareDare=p.Results.compare;

%% Do iterations

if isempty(X0)
X0=eye(size(A));
end

if isempty(S)
    S=zeros(size(Q,1),size(R,1));
end

% If S is supplied transposed
if size(S,2)==size(R,1) & size(S,1)==size(Q,1) 
S_transformed=S.';
elseif size(S,1)==size(R,1) & size(S,2)==size(Q,1) 
S_transformed=S;
end

t0=tic();
for k=0:max_iter
    
    if k==0
        Xk=X0;
    else
        Xk=Xk_plus;
    end
    
    Kk=(A*Xk*H.'+S_transformed.')/(H*Xk*H.'+R);

    Ak=A-sparse(Kk)*H;
    
    Rk=Q+A*Xk*A.'-Xk-Kk*(A*Xk*H.'+S_transformed.').'; Rk=(Rk+Rk.')/2;
    
    Nk=dlyap(Ak,Rk);
    
    Xk_plus=Xk+Nk;
    
    Xk_plus_norm=norm(Xk_plus,'fro');
    Nk_norm=norm(Nk,'fro');
    
    ratio(k+1)=Nk_norm./Xk_plus_norm;
    
    if ratio(k+1)<tol
        break
    end
    
    if k==max_iter
    disp(['***** Max iterations reached, providing solultion at k=' num2str(k) ]);
    disp(['***** Ratio Frobenius norm ||Nk|| / ||Xk_plus||  ' num2str(ratio(k+1),'%0.3e') ]);
    end
    
end
t1=toc(t0);

X_newton=Xk_plus;

%%

if DoCompareDare

t0=tic();
X_dare=idare(full(A).',full(H).',full(Q),full(R),full(S));
t2=toc(t0);

disp(['***** Calculation time Newton iterations=' num2str(t1,'%0.3f') ' s' ]);
disp(['***** Calculation time Matlab dare=' num2str(t2,'%0.3f') ' s' ]);

X_dare_norm=norm(X_dare,'fro');

X_diff_norm=norm(X_dare-X_newton,'fro');

ratio_dare=X_diff_norm./X_dare_norm;

disp(['***** Ratio Frobenius norm ||X_dare-X_newton|| / ||X_dare||=' num2str(ratio_dare,'%0.3e') ]);

end


