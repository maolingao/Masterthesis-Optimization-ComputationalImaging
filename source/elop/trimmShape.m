function mtrim = trimmShape(m,shape,fsize,offset)
msize = size(m);
switch shape
    case 'full'
        mtrim = m;
    case 'same'
        mtrim = m(ceil((fsize(1)-1)/2)+1:msize(1)-floor((fsize(1)-1)/2),ceil((fsize(2)-1)/2)+1:msize(2)-floor((fsize(2)-1)/2));
    case 'valid'
        mtrim = m((ceil((fsize(1)-1)/2))*2+1:msize(1)-(floor((fsize(1)-1)/2))*2,(ceil((fsize(2)-1)/2))*2+1:msize(2)-(floor((fsize(2)-1)/2))*2);           
    otherwise
        error('Unexpected shape in trimm.m')
end

end