function [As,Bs,Gs,Js,ys,Qs,Rs,Ss,T1,T2,T1_inv,T2_inv]=ssmod_scaleunitcov(A,B,G,J,y,Q,R,S)

%% Scale stochastic state space model so that Q=I,R=I (note: S does not become zero nor identity)
%
% Rule #1: new state xs such that x=T1*xs, i.e. state equation is multiplied by inv(T1)
% Rule #2: new output ys=inv(T2)*y, i.e. output equation is multiplied by inv(T2)
%
% Inputs:
% A: state matrix
% B: input matrix
% G: output matrix
% J: direct transmission matrix
% y: output vector
% Q: state noise covariance
% R: output noise covariance
% S: mixed noise covariance
%
% Outputs:
% As: state matrix (transformed)
% Bs: input matrix (transformed)
% Gs: output matrix (transformed)
% Js: direct transmission matrix (transformed)
% ys: output vector (transformed)
% Qs: state noise covariance (transformed)
% Rs: output noise covariance (transformed)
% Ss: mixed noise covariance (transformed)
% T1: transformation matrix, see below
% T2: transformation matrix, see below
% T1_inv: transformation matrix inverse, see below
% T2_inv: transformation matrix inverse, see below
%

%% Derivation of scaling process
%
% State equation:
% x(k+1)=A*x(k)+B*p(k)+w(k), cov(w(k))=Q;
% T1*xs(k+1)=A*T1*xs_k+B*p(k)+w(k)
% xs(k+1)=inv(T1)*A*T1*xs(k)+inv(T1)*B*p(k)+inv(T1)*w(k)
% xs(k+1)=As*xs(k)+Bs*p(k)+ws(k)
%
% where
%
% As=inv(T1)*A*T1
% Bs=inv(T1)*B
% ws(k)=inv(T1)*w(k)
%
%%%%%%%%%%%%%%%%%
% Output equation:
% y(k)=G*x(k)+J*p(k)+v(k), cov(v(k))=R;
% inv(T2)*y(k)=inv(T2)*G*T1*xs(k)+inv(T2)*J*p(k)+inv(T2)*v(k)
% ys(k)=Gs*xs(k)+Js*p(k)+vs(k)
%
% where
%
% Gs=inv(T2)*G*T1
% Js=inv(T2)*J
% vs(k)=inv(T2)*v(k)
%
%%%%%%%%%%%%%%%%%
% Next, require Qs=I
% Qs=inv(T1)*Q*inv(T1)^T
% T1*T1^T=Q
% T1*T1^T=P*Sigma*P^T
% Conclusion:
% T1=P*Sigma.^0.5
%
%%%%%%%%%%%%%%%%%
% Next, require Rs=I
% Rs=inv(T2)*R*inv(T2)^T
% T2*T2^T=R
% T2*T2^T=P*Sigma*P^T
% Conclusion:
% T2=P*Sigma.^0.5
%
%%

if ~isempty(A)
	[P1,Sigma1]=eig(Q);
	% Q=P1*S1*P1.'
	T1=(P1*Sigma1.^0.5);
	T1_inv=eye(size(T1))/T1;

	As=T1_inv*A*T1;
    
	Bs=T1_inv*B;
	Qs=T1_inv*Q*T1_inv.'; Qs=forcesym(Qs);
else 
	T1=eye(size(G,2));
	As=[];
	Bs=[];
	Qs=[];
end

if ~isempty(G)
	[P2,Sigma2]=eig(R);
	% R=P2*S2*P2.'
	T2=(P2*Sigma2.^0.5);
	T2_inv=eye(size(T2))/T2;
    
	ys=T2_inv*y;
	Gs=T2_inv*G*T1;
	Js=T2_inv*J;
	Rs=T2_inv*R*T2_inv.'; Rs=forcesym(Rs);
else 
	T2=eye(size(G,2));
	Gs=[];
	Js=[];
	Rs=[];
end

if ~isempty(A) & ~isempty(G)
	Ss=T1_inv*S*T2_inv.';
else 
	Ss=[];
end

if any(diag(Sigma1))<0; warning('Negative eigenvalues of Q'); end
if ~isreal(P1); warning('Eigenvectors of Q not real'); end

if any(diag(Sigma2))<0; warning('Negative eigenvalues of R'); end
if ~isreal(P2); warning('Eigenvectors of R not real'); end
