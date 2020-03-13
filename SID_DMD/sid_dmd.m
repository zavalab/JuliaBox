% Written by Sungho Shin (sungho.shin@wisc.edu) and Victor Zavala (victor.zavala@wisc.edu)

function [A,C,Psi,Lam]=sid_dmd(Y,n,s)

% construct data matrices in Hankel form
[m,N]=size(Y); l=N-s;
bY = zeros(m*s,l+1);;
for i=1:s
    bY((i-1)*m+1:i*m,:) = Y(:,i:l+i);
end

[P,Q] = low_rank(bY(:,1:end-1),bY(:,2:end),n); % rank-constrained regression
[A,C] = sys_mat(P,Q,m); % system matrices estimation
[Psi,Lam] = mode_dec(A,C); % mode decomposition