function [R,D] = purify_lowRank_2memlims(R1,D1,R2,D2,MEMLIM)
% low rank matrix storage limited purify, based on 2 independent purifies

Z = [R1,R2];
e = [diag(D1); diag(D2)];

[R,D] = Eig_LowRankSymmetricRealMatrix(Z,e,MEMLIM);


end