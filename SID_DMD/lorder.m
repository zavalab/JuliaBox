function [A,C]=lorder(Y,n,s,W1,W2)

if nargin==3
    W1=eye(ms)
    W2=
end

[m,N]=size(Y); N=N-s;
bYf = zeros(m*s,n); bYp= zeros(m*s,n);

for i=1:s
    W((i-1)*m+1:i*m,:) = Y(:,i:end-s-1+i);
    Wp((i-1)*m+1:i*m,:)= Y(:,i+1:end-s+i);
end

[Pt,Qt] = lrnk(W,Wp,r);
[Q,SQt,VQt] = svd(Qt,'econ');
P = Pt*VQt*SQt;

end

function [P,Q] = lrnk(bYp,bYf,n,W1,W2)
[U2,S2,V2] = svd(bYp*W2,'econ');
[U1,S1,V1] = svds(W1*bYf,n);
P=pinv(W1)*U1;
Q=U2/S2*V1*s1;
end

function [A,C] = smat(P,Q,m)
A = Q'*P;
C = P(1:m,:);
end

function [Lam,Psi] =dmode(A,C)
[Lam,Phi] = eig(A);
Psi = C*Phi;
end