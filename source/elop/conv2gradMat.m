classdef conv2gradMat < handle
    % convolution class of gradient image
    %   Detailed explanation goes here
    
    properties
        shape        % boundary condition
        GX           % convolution matrix of gradient image in x and y (u,and v) direction, in a cell array
        flag         % flag for transpose
        xsize        % size of image
        fsize        % size of blur kernel
        ysize        % size of convolved image, full case
        
        
    end
    
    methods
        % methods, including the constructor are defined in this block
        
        function obj = conv2gradMat(x,fsize,shape)
            % class constructor
            if(nargin > 0)
                if (nargin < 3)
                    obj.shape = 'full';
                else
                    obj.shape = shape;
                end
                [gx,gy]   = gradient(x);                  % gradient of x
                GXu       = conv2MatOp(gx, fsize, shape);     % convmtx of gx
                GXv       = conv2MatOp(gy, fsize, shape);     % convmtx of gy
                obj.GX    = cell(1,2);                     % cell array of convmtxs
                obj.GX{1} = GXu;
                obj.GX{2} = GXv;
                % 
                obj.xsize = GXu.xsize;
                obj.fsize = GXu.fsize;
                obj.ysize = GXu.ysize;
                obj.flag  = 0;
            end
        end
       
       function outp = mtimes(obj,f)

           switch obj.flag
               case 0 % convolution 
                   % X
                   gy       =   cellfun(@(x)x*f,  obj.GX, 'UniformOutput', false );
%                    fhat     =   obj.GX{1}'*gy{1} + obj.GX{2}'*gy{2};
                   outp     =   gy;
               case 1 % deconvolution step - correlation
                   % X^T
                   fhat     =   obj.GX{1}'*f{1} + obj.GX{2}'*f{2};
                   outp     =   fhat;
                   obj.flag =   0;
               otherwise
                   error('check class conv2gradMat');
           end
       end
       
       function obj = ctranspose(obj)
           % transpose of convolution matrix
           obj.flag = 1;
       end
       
    end
    
end

