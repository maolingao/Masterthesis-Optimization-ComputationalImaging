function [S,Y,Delta,GInv] = purify(s,y,delta,Ginv,MEMLIM)
% eliminate the component in new tuple [s,y,delta] which are
% already contained in previous spanned Krylov subspace.
% 
%
if nargin < 5
    MEMLIM = 10;
end
%     keyboard
if MEMLIM > size(s,2)       % MEMLIM > observations --> DO NOTHING
    S = s;
    Y = y;
    Delta = delta;
    GInv = Ginv;
else                        % MEMLIM < observations --> purify
    epsl = 1e-10;
    SIGMA = 'LM';           
    G = s'*y;
    [U,D] = eigs(G,(MEMLIM),SIGMA);     % preserve G's eigenvalue with MEMLIM largest magnetitude 
                                        % and cooresponding eigenvector
    %
    % [U,D] = eigs(Ginv,(MEMLIM),SIGMA);
    %
    % keyboard
    U = real(U(:,1:MEMLIM));
    D = real(D(1:MEMLIM,1:MEMLIM));
    %
%     display(sprintf('Dmin=--------------------------------------%d',min((diag(D)))));
%     U = bsxfun(@times,U,((sqrt(diag(D))+eps).\1)');    % make D unit matrix, numerical stabil if invert it
%     D = eye(MEMLIM);
%
%--------------- unittest_1 ---------------%
% should give a matrix U'*U approximately = I, diagnal matrix(actually inv(D))
%     keyboard
%     figure(3), imagesc(real((U'*U))), colormap gray, axis image, colorbar('southoutside')
%---------------    END     ---------------%
%
    %
    S = s*U;                            % new tuple of {S, Y, Delta}
    Y = y*U;                            % S should be orthogonal to Y, S'*Y = D
    Delta = delta * U;
%
%--------------- unittest_2 ---------------%
% should give a matrix with amost all zero element
%     keyboard
%     S'*Y - D
%---------------    END     ---------------%
%
    GInv = diag((diag(D) + eps).\1);
end
end
%
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