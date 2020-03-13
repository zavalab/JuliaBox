% Written by Sungho Shin (sungho.shin@wisc.edu) and Victor Zavala (victor.zavala@wisc.edu)

function [P,Q] = low_rank(bYp,bYf,n)

[U2,S2,V2]=svd(bYp,'econ');
[P,S1,V1]=svds(bYf*V2,n);
Q = U2/S2*V1*S1;