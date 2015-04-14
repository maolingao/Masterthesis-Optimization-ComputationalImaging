function [R,D] = purify_lowRank_obs(s,y,delta,MEMLIM)
% low rank matrix storage limited purify, based on observations within
% linear problems
G = s' * y;
G = 1/2*(G + G');

Ginv = pinv(G);
V = s * Ginv;
% V = s / G;
U = delta - 1/2 * V * y' * delta;

Z = 1/sqrt(2) * [V - U, V + U];
i = ones(size(V,2),1);
e = [-i; i];
% E = diag(-i,i); Ztilde = E * Z';

[R,D] = Eig_LowRankSymmetricRealMatrix(Z,e,MEMLIM);



end