% Written by Sungho Shin (sungho.shin@wisc.edu) and Victor Zavala (victor.zavala@wisc.edu)

function [A,C] = sys_mat(P,Q,m)
A = Q'*P;
C = P(1:m,:);
