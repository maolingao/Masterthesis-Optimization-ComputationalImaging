classdef conv2MatOp < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x           % matrix, image
        xsize       % size of image
        f           % matrix blur kernel
        fsize       % size of blur kernel
        ysize       % size of convolved image, full case
        shapeSize   % size of convolved image according shape
        yfull       % full convolved image
        shape       % shape option
        flag=0;     % flag for devonvolution
        offset;     % offset for padding
        fft_x;      % cache of x's fft
        Px          % pad and trim mtx for x
        Pf          % pad and trim mtx for f
        Py          % pad and trim mtx for y
        
    end
    
    methods
        % methods, including the constructor are defined in this block
        
        function obj = conv2MatOp(x,fsize,shape)
            % class constructor
            if(nargin > 0)                
                obj.x = x;
                obj.xsize = size(x);
                obj.fsize = fsize;
                obj.ysize = obj.xsize + obj.fsize - ones(1,2);
                if (nargin < 3)
                    obj.shape = 'full';
                else
                    obj.shape = shape;
                end
                switch shape
                    case 'full'
                        obj.shapeSize = obj.ysize;
                    case 'same'
                        obj.shapeSize = max(obj.xsize,obj.fsize);
                    case 'valid'
                        obj.shapeSize = max(obj.xsize,obj.fsize) - min(obj.xsize,obj.fsize) + 1;
                    otherwise
                        error('check conv2MatOp shape')
                end
                obj.offset = [0 0];
                obj.Px = patimat(obj.ysize,obj.xsize,obj.offset);
                obj.Pf = patimat(obj.ysize,obj.fsize,obj.offset);
                kernelSize = min(obj.xsize, obj.fsize); % in trimmShape.m, it is actually fsize 
                obj.Py = patimat(obj.shape,obj.ysize,kernelSize,obj.offset);
                obj.fft_x = fft2(obj.Px*x);
                
%                 obj = buildConvMtx(obj);
            end
        end
       
       function outp = mtimes(obj,f)

           switch obj.flag
               case 0 % convolution 
                   obj.yfull = real(ifft2(obj.fft_x.*fft2(obj.Pf*f)));  
%                    figure, imagesc(obj.yfull)
                   outp = obj.Py'*obj.yfull;
%                    figure, imagesc(outp)

               case 1 % deconvolution
%                    figure,imagesc(padm(obj,f,obj.ysize,offset_centre))
                   yShape = f;
                   ffull = real(ifft2(conj(obj.fft_x).*fft2(obj.Py*yShape))); % F transpose                 
%                    figure, imagesc(ffull)
                   outp = obj.Pf'*ffull;
%                    figure, imagesc(outp)
                   obj.flag = 0;
               otherwise
                   error('check class conv2MatOp');
           end   
%            outp = obj.yfull;
       end
       
       function obj = ctranspose(obj)
           % transpose of convolution matrix
           obj.flag = 1;
       end
       
    end
    
end

