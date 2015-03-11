function [R,D] = driftH(R1,D1,R2,D2,MEMLIM,alpha)
% calculate H_3_prior = H_2 + alpha * (H_2 - H_1)
% H_3_prior = R * D * R'
% scaled drift of H from two linear problems

alpha = alpha^(1/3);
Z = [R2,       alpha*R2,        alpha*R1];
e = [diag(D2); alpha*diag(D2); -alpha*diag(D1); ];

[R,D] = Eig_LowRankSymmetricRealMatrix(Z,e,MEMLIM);




end