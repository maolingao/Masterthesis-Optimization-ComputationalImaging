classdef hessianMatrix < handle
% infered hessian matrix class for probabilistic solver
% Detailed explanation goes here
properties
        H           % infered hessian matrix
        s           % corrected step error
        y           % corrected step residual
        delta       % difference btw corrected step error and back-calculated corrected step residual
        i           % iteration index
        Wfun        % W prior covariance function
        H0fun       % H prior mean function
        R           % dominant eigenvectors
        D           % dominant eigenvalues
end
methods
    function obj = hessianMatrix(H,s,y,delta,R,D,Wfun,H0fun)
        if nargin > 0
            obj.H = H;
            if exist('s','var') && ~isempty(s)
                obj.s = s;
                obj.i = size(s,2) + 1;
            else
                obj.i = 1;
            end
            if exist('y','var') && ~isempty(y)
                obj.y = y;
            end
            if exist('delta','var') && ~isempty(delta)
                obj.delta = delta;
            end
            if exist('R','var') && ~isempty(R)
                obj.R = R;
            end
            if exist('D','var') && ~isempty(D)
                obj.D = D;
            end
            if exist('Wfun','var') && ~isempty(Wfun)
                obj.Wfun = Wfun;
            elseif exist('R','var') && ~isempty(R)
                obj.Wfun = @(x) R * (D * (R' * x));
            else
                obj.Wfun = @(x) x;
            end
            if exist('H0fun','var') && ~isempty(H0fun)
                obj.H0fun = H0fun;
            elseif exist('R','var') && ~isempty(R)
                obj.H0fun = @(x) R * (D * (R' * x));
            else
                obj.H0fun = @(x) x;
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
        xresp = reshape(x,[size(obj.H,1),numel(x)./size(obj.H,1)]);
        
        % the part of current observations
        if ~isempty(obj.s)
           % compute W*y via handle
            if isa(obj.Wfun,'function_handle')
                z = obj.Wfun(obj.y);
            else
                error('malformed covariance function')
            end
            G = obj.y' * z;
            %----------------------------%
            % pseudo-inverse
            Ginv = pinv(G);
            ZX = (z' * xresp);
            GinvZX = Ginv * ZX;
            ZGinv = z * Ginv;
            %----------------------------%
            % backslash
%             ZGinv   =  (G \ z')';
%             GinvZX  =  ZGinv' * xresp;
            %----------------------------%
            tail = obj.delta * GinvZX + ZGinv * (obj.delta' * xresp) - (ZGinv * (obj.y' * obj.delta)) * GinvZX;
        else
            tail = 0;
        end
        % *** update form: H = H0 + tail ***
        outp = obj.H0fun(xresp) + tail;   % output = vec(I*x) + tail
    end
    
    function outp = times(obj,x)
    % multiplication with matrix and vector(last term of H)
        epsl = 1e-30;
        tail = 0; % (update)
        % only use the last term of s, y and Delta to update H
        den = (obj.s(:,end)'*obj.y(:,end) );
        tail = tail + (den + epsl )\(obj.delta(:,end)*(obj.s(:,end)'*vec(x)) ) + ...
        (den + epsl )\(obj.s(:,end)*(obj.delta(:,end)'*vec(x)) ) - ...
        (den^2 + epsl )\(obj.s(:,end)*(obj.delta(:,end)'*obj.y(:,end))*(obj.s(:,end)'*vec(x)) );
        outp = tail; % output = tail(last term)
    end
end
end