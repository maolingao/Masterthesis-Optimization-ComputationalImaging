function [R,D] = purify_lowRank_2memlims(R1,D1,s,y,delta,MEMLIM,option)
% low rank matrix storage limited purify, based on information from
% previous accumulated knowledge and new observations.

[V,U] = rank2form(s,y,delta,option);
r = 1/sqrt(2) * [V - U, V + U];
i = ones(size(V,2),1);
d = [-i; i];


Z = [R1,r];
e = [diag(D1); d];

[R,D] = Eig_LowRankSymmetricRealMatrix(Z,e,MEMLIM);


end