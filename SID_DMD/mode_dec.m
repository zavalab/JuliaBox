% Written by Sungho Shin (sungho.shin@wisc.edu) and Victor Zavala (victor.zavala@wisc.edu)

function [smode,tmode] = mode_dec(A,C)

[Phi,Lam] = eig(A);
smode = C*Phi;
tmode = diag(Lam);