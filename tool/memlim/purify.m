function [S,Y,Delta,GInv] = purify(s,y,delta,MEMLIM,lambda)
% eliminate the component in new tuple [s,y,delta] which are
% already contained in previous spanned Krylov subspace.
%
%
if nargin < 4
MEMLIM = 10;
end
if nargin < 5
lambda = 0;
end

if MEMLIM > size(s,2) % MEMLIM > observations --> guarentee s and y are conjugate by eig-decomp
    G       =   s'*y;
    [U,D]   =   eig(G);
    idx     =   sum(diag(D)>1e-6);
    U       =   real(U(:,1:idx));
    D       =   lambda * clip(real(D(1:idx,1:idx)),inf,0);
    S       =   s * U * sqrt(lambda);
    Y       =   y * U * sqrt(lambda);
    Delta   =   delta * U * sqrt(lambda);
    GInv    =   diag((diag(D) + eps).\1);
else % MEMLIM < observations --> purify
    epsl    =   1e-10;
    SIGMA   =   'LM';
    G       =   s'*y;

    [U,D]   =   eigs(G,(MEMLIM),SIGMA); % preserve G's eigenvalue with MEMLIM largest magnetitude
    % and cooresponding eigenvector
   
    idx     =   sum(diag(D)>1e-6);
    if idx == 0
        S       =   [];
        Y       =   [];
        Delta   =   [];
        return
    end
    if idx >= MEMLIM
        U       =   real(U(:,1:MEMLIM));
        D       =   clip(real(D(1:MEMLIM,1:MEMLIM)),inf,0);
    else
        U       =   real(U(:,1:idx));
        D       =   clip(real(D(1:idx,1:idx)),inf,0);
    end

    %--------------- unittest_1 ---------------%
    % should give a matrix U'*U approximately = I, diagnal matrix(actually inv(D))
    % keyboard
    % figure(3), imagesc(real((U'*U))), colormap gray, axis image, colorbar('southoutside')
    %--------------- END ---------------%
    %
    D       =   lambda * D; % memory strength
    S       =   s * U * sqrt(lambda); % new tuple of {S, Y, Delta}
    Y       =   y * U * sqrt(lambda); % S should be orthogonal to Y, S'*Y = D
    Delta   =   delta * U * sqrt(lambda);

    %--------------- unittest_2 ---------------%
    % should give a matrix with amost all zero element
    % keyboard
    % S'*Y - D
    %--------------- END ---------------%
    %
    GInv    =   diag((diag(D) + eps).\1);
end
end
%
% unittest
%{
n = 80;
u = rand(n,1);
Q = RandomRotation(n);
D = diag(u);
A = Q*D*Q'; % symmetric pos def matrix with eigenvalue btw [0,1]
[V,D] = eig(A); % V eigenvector, mutally orthonormal
%
s = bsxfun(@times,V(:,1:30),rand(1,30));
G = s'*s; % symmetric diagnal matrix
Ginv = G\eye(size(G));
delta = zeros(size(s));
MEMLIM = 10;
lambda = 0;
[S,Y,Delta,GInv] = purify(s,s,delta,Ginv,MEMLIM,lambda);
%
%}
