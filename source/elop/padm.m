function [mpad] = padm(m,fullsize,offset) 
% pad-matrix for a matrix m


mpad = zeros(fullsize);
msize = size(m);
%         fsizehat = fullsize - msize;
mpad(1:msize(1),1:msize(2)) = m;        
%         mpad = circshift(mpad,ceil(fsizehat./2));
        
end