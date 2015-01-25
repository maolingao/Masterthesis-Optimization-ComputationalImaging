function [S,Y,Delta,GInv] = purify(s,y,delta,Ginv,MEMLIM)
% eliminate the component in new tuple [s,y,delta] which are
% already contained in previous spanned Krylov subspace.
% 
%
if nargin < 5
    MEMLIM = 10;
end
epsl = 1e-10;
SIGMA = 'LM';
[U,D] = eigs(Ginv,MEMLIM,SIGMA);
U = bsxfun(@rdivide, U, max(abs(U)));
S = s*U;
Y = y*U;
keyboard
U'*U;
Delta = delta * U;
GInv = (S'*Y + epsl)\eye(size(S,2));
end

% unittest
%{
n = 80;
u = rand(n,1);
Q = RandomRotation(n); 
D = diag(u);
A = Q*D*Q';         % symmetric pos def matrix with eigenvalue btw [0,1]
[V,D] = eig(A);     % V eigenvector, mutally orthonormal
%
s = bsxfun(@times,V(:,1:30),rand(1,30));
G = s'*s;           % symmetric diagnal matrix 
Ginv = G\eye(size(G));
delta = zeros(size(s));
[S,Y,Delta,GInv] = purify(s,s,delta,Ginv);
%
%}