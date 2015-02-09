classdef hessianMatrix < handle
    % infered hessian matrix class for probabilistic solver
    %   Detailed explanation goes here
    
    properties
        H           % infered hessian matrix
        s           % corrected step error
        y           % corrected step residual
        delta       % difference btw corrected step error and back-calculated corrected step residual
        i           % iteration index
        Ginv0       % inverse of last Gramm matrix
        scale       % scale of diagonal element in H
    end
    
    methods
        function obj = hessianMatrix(H,s,y,delta,Ginv0,i)
            if nargin > 0
                obj.H = H;
                obj.scale = unique(diag(obj.H));
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
                if exist('Ginv0','var')
                    obj.Ginv0 = Ginv0;
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
%                 keyboard
                % working version, using all terms to update H
                %{
                tail = 0; % (update)
                for k = 1 : obj.i - 1
                    den = (obj.s(:,k)'*obj.y(:,k) );
                    tail = tail + (den + epsl )\(obj.delta(:,k)*(obj.s(:,k)'*vec(x))  ) + ...
                        (den + epsl )\(obj.s(:,k)*(obj.delta(:,k)'*vec(x))  ) - ...
                        (den^2 + epsl )\(obj.s(:,k)*(obj.delta(:,k)'*obj.y(:,k))*(obj.s(:,k)'*vec(x)) );
                end
                %}
                % debug, using full H
                %{%
                if ~isempty(obj.s)
                    if size(obj.s,2) == size(obj.Ginv0,1)
                        Ginv = obj.Ginv0;
                    elseif size(obj.s,2) < 2 
%                         keyboard
                        Ginv = 1/((obj.s)'*obj.y + epsl);
                    else
%                         keyboard
                        Ginv = invGram(obj.Ginv0,obj.s,obj.y);
                    end
%                     SGinv = obj.s * Ginv;
%                     SGinvDelta = SGinv * (obj.delta)';
%                     tailM = SGinvDelta + SGinvDelta' - SGinv * obj.y' * SGinvDelta';
%                     tail = tailM * vec(x);                    
                    %----------------------------%
%                     keyboard
                    SX = (obj.s' * vec(x));
                    GinvSX = Ginv * SX;
                    SGinv = obj.s * Ginv;
                    tail = obj.delta * GinvSX + SGinv * (obj.delta' * vec(x)) - SGinv * (obj.delta' * obj.y) * GinvSX;
                    %----------------------------%
                    obj.Ginv0 = Ginv;
                else
                    tail = 0;
                end
                %}                
                %{
                if ~isempty(obj.s)
                    if size(obj.s,2) == size(obj.Ginv0,1)
%                         keyboard
                        Ginv = obj.Ginv0;
                    elseif size(obj.s,2) < 2 
                        Ginv = 1/((obj.s)'*obj.s + epsl);
                    else
%                         keyboard
                        M = (obj.s)'*obj.y;
                        Ginv = (M'*M) \ eye(size(obj.s,2)) * M';
                    end
                    SGinv = obj.s * Ginv;
                    SGinvDelta = SGinv * (obj.delta)';
                    tailM = SGinvDelta + SGinvDelta' - SGinv * obj.y' * SGinvDelta';
                    tail = tailM * vec(x);
                    obj.Ginv0 = Ginv;
                else
                    tail = 0;
                end
                %}
                % *** update form: H = H0 + tail ***
                outp = obj.scale * vec(x) + tail;              % output = vec(I*x) + tail
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