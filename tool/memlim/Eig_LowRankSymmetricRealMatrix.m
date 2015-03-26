function [R,D] = Eig_LowRankSymmetricRealMatrix(U,e,cutoff)
% [R,D] = Eig_LowRankSymmetricRealMatrix(U,e,cutoff)
% returns the 'economical' eigen-decomposition of the symmetric real-valued
% matrix of low rank     A = U * E * U'       where E = diag(e);
% That is, A = R * D * R' with a rectangular matrix R such that R' * R = Id
% Only eigenvectors corresponding to nonzero-eigenvalues are returned.
%
% If E is a vector, it is treated as the elements of a diagonal matrix.
%
% IF cutoff < 1, only the
% eigen-dimensions corresponding to eigenvalues lambda_i with
% lambda_i > cutoff * max(lambda) are returned (note that there is no
% absolute value here. So this can be used to truncate to positive definite
% matrices. 
% 
% IF cutoff >= 1 and an integer, cutoff=M, only the M largest eigenvalues are returned.
%
% The output is sorted in descending order by the _absolute_ value of the eigenvalues.
%
% 
%
% (c) 2014, Philipp Hennig
% 2015 modified, Maolin Gao

if any(~isreal(U(:))) || any(~isreal(e(:)))
    error ('this operation only works for real-valued matrices');
end
% assert(size(U,1) > size(U,2)); assert(isvector(e))
if (size(U,1) > size(U,2))
    if isrow(e); e = e'; end

    [V,D] = eig(bsxfun(@times,e,(U' * U)));   % eigendecomposition of complement

    % if any(~isreal(diag(D)))
    %     warning('Found complex eigenvalues of a symmetric real matrix. This indicates the operation is not stable. Imaginary parts will be dropped at end of calculation.')
    % end

    R     = U * V * D;                        % unnormalized eigenvectors of A
    d     = diag(D);
else
    warning('[Eig_LowRankSymmetricRealMatrix.m] : more columns than rows, maybe computational costly.')
    [R,D] = eig(U * diag(e) * U');   % eigendecomposition of complement
    d     = diag(D);
end

R = real(R); d = real(d);                    % this can help stabilize. 
R = bsxfun(@rdivide,R,sqrt(sum(R.^2,1)));    % normalize eigenvectors

[~,seq] = sort(d,'descend');     % sort
d       = d(seq);
R       = R(:,seq);              % make sure R has right order

idx = d > 1e-6;
d = d(idx);
R = R(:,idx);

if nargin > 2          % remove small eigenvalues below cutoff    
    if cutoff < 1
        R       = R(:,d > cutoff * max(d));
        d       = d(d > cutoff * max(d));
    elseif round(cutoff) - cutoff == 0
        cutoff = min([cutoff, length(d)]);
        R       = R(:,1:cutoff);
        d       = d(1:cutoff);
    else error('malformed input cutoff. Must be float < 1, or integer > 1.')
    end
end
D       = diag(d);

    
end

%%% unit test:
% N = 10; M = 5;
% U = randn(N,M); e = randn(M,1);  A = U * diag(e) * U';
% [R,D] = Eig_LowRankSymmetricRealMatrix(U,e);
% A * R - R * D
% ^^^^ should give a NxM matrix of almost zeros.

