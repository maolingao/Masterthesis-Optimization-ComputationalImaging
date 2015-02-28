function [b,c] = factorize(a)
% factorize scalar a intp product of b and c : a = b*c
% which b and c and most close to each other, and b < c

if ~exist('a','var') || a < 1
    display('[factorize.m] : number to factorize should be larger than 1.')
    return
else
    b       =   [];
    f       =   factor(a);              % prime factors
    aroot   =   sqrt(a);
    if round(aroot) == aroot, b = aroot; c = aroot; return; end  % if square root is integer, return b=c=sqrt(a);
    for i   =   1 : length(f)
        C   =   nchoosek(f,i);          % Combinatorics, Binomial coefficient or all combinations
        C   =   unique(C,'rows');       % unique combinations
        b   =   [b; prod(C, 2)];        % possible products of all unique combinations
    end
    [~,ind] = min(abs(b-aroot));        % b value closest to sqrt(a);
    b       =   b(ind);
    c       =   a/b;    
    if b < c                            % ensure output oder : b < c
        return
    else
        [b,c] = swap(b,c);
    end
end

end