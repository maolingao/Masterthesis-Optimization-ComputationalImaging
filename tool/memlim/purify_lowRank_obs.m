function [R,D] = purify_lowRank_obs(s,y,delta,MEMLIM,option)
% low rank matrix storage limited purify, based on observations within
% linear problems

[V,U] = rank2form(s,y,delta,option);

r = 1/sqrt(2) * [V - U, V + U];
i = ones(size(V,2),1);
d = [-i; i];
% E = diag(-i,i); Ztilde = E * Z';

[R,D] = Eig_LowRankSymmetricRealMatrix(r,d,MEMLIM);



end