function [mpad] = padmShape(m,shape,fullsize,offset) 
% pad-matrix for a matrix m

switch shape
    case 'full'
        mpad = m;
    case 'same'
        mpad = zeros(fullsize);
        msize = size(m);
        fsizehat = fullsize - msize;
        mpad(1:msize(1),1:msize(2)) = m;        
        mpad = circshift(mpad,ceil(fsizehat./2));
    case 'valid'
        mpad = zeros(fullsize);
        msize = size(m);
        fsizehat = fullsize - msize;
        mpad(1:msize(1),1:msize(2)) = m;        
        mpad = circshift(mpad,ceil(fsizehat./2));
    otherwise
        error('Unexpected shape in padmShape.m');
        
end

end