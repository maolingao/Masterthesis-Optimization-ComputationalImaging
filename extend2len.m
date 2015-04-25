function x_extended = extend2len(x,len)

if ~exist('len','var'); x_extended = x; return; end
if size(x,1) < size(x,2); x = x'; end

lx = size(x,1);
lenDiff = len - lx;
if lenDiff == 0; x_extended = x; return; end

if lenDiff > 0
    x_append = x(end)*ones(lenDiff,1);
    x_extended = [x;x_append];
else
    error('in [extend2len.m] : len should be larger than signal.')
end
end