function mhat = trimm(m,fsize,offset)
% crop m to recover original m size

% msize = size(m);
% mhatsize = msize - fsize + ones(1,2);
% mhat = m(1:mhatsize(1),1:mhatsize(2));
mhat = m(1:fsize(1),1:fsize(2));

end