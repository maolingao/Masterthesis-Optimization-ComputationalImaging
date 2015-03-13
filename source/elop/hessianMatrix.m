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
        R           % dominant eigenvectors
        D           % dominant eigenvalues
    end
    
    methods
        function obj = hessianMatrix(H,s,y,delta,R,D,Ginv0)
            if nargin > 0
                obj.H = H;
                if exist('s','var') && ~isempty(s)
                    obj.s = s;
                    obj.i = size(s,2) + 1;
                else
                    obj.i = 1;
                end
                if exist('y','var')  && ~isempty(y)
                    obj.y = y;
                end
                if exist('delta','var') && ~isempty(delta)
                    obj.delta = delta;
                end   
                if exist('R','var')  && ~isempty(R)
                    obj.R = R;
                end
                if exist('D','var')  && ~isempty(D)
                    obj.D = D;
                end
                if exist('Ginv0','var') && ~isempty(Ginv0)
                    obj.Ginv0 = Ginv0;
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
                % debug, using full H, pseudo-invers
                % the part of previous observations, storage limited
                if ~isempty(obj.R) && ~isempty(obj.D)
                    tail = obj.R * (obj.D * (obj.R' * vec(x)));
                else 
                    tail = 0;
                end
                % the part of current observations
                if ~isempty(obj.s)
%                     G = obj.s'*obj.y;
                    G = obj.y' * obj.y;
                    G = 1/2 * ( G + G');
                    Ginv = pinv(G);
                    %----------------------------%
                    SX = (obj.s' * vec(x));
                    GinvSX = Ginv * SX;
                    SGinv = obj.s * Ginv;
%                     SGinv = obj.s / G;
%                     GinvSX = SGinv' * vec(x);
                    tail = tail + obj.delta * GinvSX + SGinv * (obj.delta' * vec(x)) - SGinv * (obj.delta' * obj.y) * GinvSX;
                    %----------------------------%
                    obj.Ginv0 = Ginv;
                else
                    tail = tail + 0;
                end
                % *** update form: H = H0 + tail ***
                outp = vec(x) + tail;              % output = vec(I*x) + tail
                % ########################
                
        end
        
        function outp = times(obj,x)
            % multiplication with matrix and vector(last term of H)
                epsl = 1e-30;
                tail = 0; % (update)
                % only use the last term of s, y and Delta to update H
                % aim : speed up
                    den = (obj.s(:,end)'*obj.y(:,end) );
                    tail = tail + (den + epsl )\(obj.delta(:,end)*(obj.s(:,end)'*vec(x)) ) + ...
                        (den + epsl )\(obj.s(:,end)*(obj.delta(:,end)'*vec(x)) ) - ...
                        (den^2 + epsl )\(obj.s(:,end)*(obj.delta(:,end)'*obj.y(:,end))*(obj.s(:,end)'*vec(x))  );

                outp = tail; % output = tail(last term)
                
        end
    end
end