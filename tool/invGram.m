function Ginv = invGram(Ginv0,S,Y)
% invert Gramm matrix G = S'*Y
%
epsl = 1e-9;
s = S(:,end);               % last observation in S
y = Y(:,end);               % last observation in Y
a = S(:,1:end-1)'*y;        % new column in G apart from the last element
c = s'*y;                   % last element of new column in G
%
keyboard
M = (c - a'*Ginv0*a + epsl)\1;     % last new element in Ginv
Greg = a' * Ginv0;          % most expensive computation
lowerLeft = -M * Greg;      % lower-left block in Ginv
%
Ginv = [Ginv0 + Ginv0 * a * M * Greg , lowerLeft';
        lowerLeft                    , M         ];


end

% unittest
% show give a matrix with amost zero elements
%{
n = 80;
u = rand(n,1);
Q = RandomRotation(n); 
D = diag(u);
A = Q*D*Q';
[V,D] = eig(A);
%
S = bsxfun(@times,V(:,1:30),rand(1,30));
G0 = S'*S;
Ginv0 = G0\eye(size(G0));
%
S = [S,V(:,31)*2];
G = S'*S;
Ginv = invGram(Ginv0,S,S);
Ginv - G\eye(size(G))
%}