classdef hessianMatrix < handle
    % infered hessian matrix class for probabilistic solver
    %   Detailed explanation goes here
    
    properties
        H           % infered hessian matrix
        s           % corrected step error
        y           % corrected step residual
        delta       % difference btw corrected step error and back-calculated corrected step residual
        i           % iteration index
        l           % low rank update matrix H = H0 + l*r
        r           % low rank update matrix H = H0 + l*r
        
    end
    
    methods
        function obj = hessianMatrix(H,l,r,s,y,delta,i)
            if nargin > 0
                obj.H = H;
                if exist('l','var')
                    obj.l = l;
                end
                if exist('r','var')
                    obj.r = r;
                end
                if exist('s','var')
                    obj.s = s;
                    obj.i = size(s,2) + 1;
                end
                if exist('y','var')
                    obj.y = y;
                end
                if exist('delta','var')
                    obj.delta = delta;
                end   
                if exist('i','var')
                    obj.i = i;
                else
                    obj.i = 1;
                end             
             end
            
        end
        
        function obj = plus(obj,s,y,delta)
            % update H je iteration
            obj.s(:,obj.i) = vec(s);
            obj.y(:,obj.i) = vec(y);
            obj.delta(:,obj.i) = vec(delta);
            obj.i = obj.i + 1;
        end
        
        function outp = mtimes(obj,x)
            % multiplication with matrix and vector
            
                epsl = 1e-30;
                X = x;
%                 keyboard
                for col_ind = 1:size(X,2)
                    x = X(:,col_ind);
                % working version, using all terms to update H
                %{%
                tail = 0; % (update)
                for k = 1 : obj.i - 1
%                 tail = tail + (obj.s(:,k)'*obj.y(:,k) + epsl )\(obj.delta(:,k)*(obj.s(:,k)'*vec(x)) + epsl ) + ...
%                     (obj.s(:,k)'*obj.y(:,k) + epsl )\(obj.s(:,k)*(obj.delta(:,k)'*vec(x)) + epsl ) - ...
%                     ((obj.s(:,k)'*obj.y(:,k))^2 + epsl )\(obj.s(:,k)*(obj.delta(:,k)'*obj.y(:,k))*(obj.s(:,k)'*vec(x)) + epsl );
%                 end
                    den = (obj.s(:,k)'*obj.y(:,k) );
                    tail = tail + (den + epsl )\(obj.delta(:,k)*(obj.s(:,k)'*vec(x))  ) + ...
                        (den + epsl )\(obj.s(:,k)*(obj.delta(:,k)'*vec(x))  ) - ...
                        (den^2 + epsl )\(obj.s(:,k)*(obj.delta(:,k)'*obj.y(:,k))*(obj.s(:,k)'*vec(x)) );
                end
                %}
                % ########################
                % check, update fasion of H
                %{%
                % *** update form: H = H0 + l*r' + tail ***
                if ~isempty(obj.l)
                    lrt = obj.l*(obj.r'*vec(x));
                    outp(:,col_ind) = vec(x) + vec(lrt) + tail; % output = vec(I*x) + l*r'*vex(x) + tail
                else
                    outp(:,col_ind) = vec(x) + tail;            % output = vec(I*x) + tail
                end
                end
                %}
                % *** update form: H = H0 + tail ***
%                 outp = vec(x) + tail;              % output = vec(I*x) + tail
                % ########################
                
        end
        
        function outp = times(obj,x)
            % multiplication with matrix and vector(last term of H)
                epsl = 1e-30;
                tail = 0; % (update)
                % only use the last term of s, y and Delta to update H
                % aim : speed up
                % prob : residual after step from Ax=b to Ax'=b' not any more conjugate with all old search directions
%                 endInd = obj.i - 1;
%                 if endInd > 0
%                 tail = tail + (obj.s(:,end)'*obj.y(:,end) + epsl )\(obj.delta(:,end)*(obj.s(:,end)'*vec(x)) + epsl ) + ...
%                     (obj.s(:,end)'*obj.y(:,end) + epsl )\(obj.s(:,end)*(obj.delta(:,end)'*vec(x)) + epsl ) - ...
%                     ((obj.s(:,end)'*obj.y(:,end))^2 + epsl )\(obj.s(:,end)*(obj.delta(:,end)'*obj.y(:,end))*(obj.s(:,end)'*vec(x)) + epsl );

                    den = (obj.s(:,end)'*obj.y(:,end) );
                    tail = tail + (den + epsl )\(obj.delta(:,end)*(obj.s(:,end)'*vec(x)) ) + ...
                        (den + epsl )\(obj.s(:,end)*(obj.delta(:,end)'*vec(x)) ) - ...
                        (den^2 + epsl )\(obj.s(:,end)*(obj.delta(:,end)'*obj.y(:,end))*(obj.s(:,end)'*vec(x))  );
%                 end

                outp = tail; % output = tail(last term)
%                 outp = vec(obj.H*x) + tail; % H_i+1 <-- H_i + (update)
                
        end
    end
end