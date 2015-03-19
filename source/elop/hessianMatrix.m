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
            end
            if exist('H0fun','var') && ~isempty(H0fun)
                obj.H0fun = H0fun;
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
        % the part of previous observations, storage limited
        if ~isempty(obj.R) && ~isempty(obj.D)
            tail = obj.R * (obj.D * (obj.R' * vec(x)));
        else
            tail = 0;
        end
        % the part of current observations
        if ~isempty(obj.s)
           % compute W*y via handle
            if isa(obj.Wfun,'function_handle')
                z = obj.Wfun(obj.y);
            elseif isa(obj.Wfun,'char') && strcmp(obj.Wfun,'BFGS')
            % implicit W=H function:
                z = obj.s;
            elseif isa(obj.Wfun,'char') && strcmp(obj.Wfun,'GS')
            % explicitly W=H0 function:
                if ~isempty(obj.R) && ~isempty(obj.D)
                    z = obj.R * obj.D * (obj.R' * obj.y);
                else
                    z = obj.y;
                end
            else
                error('malformed covariance function')
            end
            G = obj.y' * z;
            %----------------------------%
            % pseudo-inverse
            Ginv = pinv(G);
            ZX = (z' * vec(x));
            GinvZX = Ginv * ZX;
            ZGinv = z * Ginv;
            %----------------------------%
            % backslash
%             ZGinv   =  (G \ z')';
%             GinvZX  =  ZGinv' * vec(x);
            %----------------------------%
            tail = tail + obj.delta * GinvZX + ZGinv * (obj.delta' * vec(x)) - (ZGinv * (obj.y' * obj.delta)) * GinvZX;
        else
            NOP;
        end
        % *** update form: H = H0 + tail ***
        outp = obj.H0fun(vec(x)) + tail;   % output = vec(I*x) + tail
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