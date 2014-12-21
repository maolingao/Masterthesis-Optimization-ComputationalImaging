classdef patimat < handle
    % patimat Summary of this class goes here
    %   Detailed explanation goes here
    % pad(pa) and trim(ti) matrix(mat) operator
    
    properties
        shape               % conv boundary condition
        fullsize            % full size
        targetSize          % size want to have, at trimm.m / in trimmShape.m, it is actually fsize
        offset              % offset
        padFlag             % flag for choosing pad operation: padm.m/padmShape.m
        trimFlag            % flag for trim.m
    end
    
    methods
        % methods, including the constructor are defined in this block
        function obj = patimat(varargin)    % (shape,fullsize,targetSize,offset)
            % class constructor
            if nargin == 4
                obj.shape = varargin{1};
                obj.fullsize = varargin{2};
                obj.targetSize = varargin{3};
                obj.offset = varargin{4};  
                obj.padFlag = 1;        % pad r.t. shape
            elseif nargin == 3
                obj.fullsize = varargin{1};
                obj.targetSize = varargin{2};
                obj.offset = varargin{3};   
                obj.padFlag = 0;        % pad no matter shape
            end   
                obj.trimFlag = 0;
                
            
        end
        %mtimes
        function outp = mtimes(obj,m)
            switch obj.trimFlag
                case 0
                    if obj.padFlag == 0
                        mpad = padm(m,obj.fullsize,obj.offset);
                    elseif obj.padFlag == 1
                        mpad = padmShape(m,obj.shape,obj.fullsize,obj.offset);
                    else
                        error('Unexpected padFlag in patimat')
                    end
                    outp = mpad;
                case 1
                    if obj.padFlag == 0
                        mhat = trimm(m,obj.targetSize,obj.offset);
                        outp = mhat;
                    elseif obj.padFlag == 1
                        mtrim = trimmShape(m,obj.shape,obj.targetSize,obj.offset);
                        outp = mtrim;
                    else
                        error('Unexpected padFlag in patimat')
                    end
                    obj.trimFlag = 0;
                otherwise
                    error('check class patimat')
            end
        end
        % ctranspose
        function obj = ctranspose(obj)            
           obj.trimFlag = 1;
        end
    end
    
end

