function [x,H,x_stack,residual] = pncg_Hmfd(A,b,H,x_start,tol,iter,option)
% probabilistic solver, conjugate gradient
% written for mfd_test script, to verify idea of multi-frame deconv(update of H with class hessianMatrix)

% example - input x_start = x, output must stay at x_start
startup;

%{

figure(999), clf
figure(997), clf
n = 60;
u =-10 * log(rand(n,1));u(1:5) = 100*u(1:5);
Q = RandomRotation(n); 
D = diag(u);
A = Q*D*Q';
% figure(5), imagesc(A),colormap
x = rand(n,1);
b = A*x;
H  = eye(size(A));
x_start = x;
tic, [x,H] = pncg_Hmfd(A,b,H,x_start); toc
%}
if nargin < 3
    H = hessianMatrix(eye(size(A)));
end
if nargin < 4
    x_start = b;
end
if nargin < 5
    tol = 10^-10;
end
if nargin < 6
    iter = 100;
end
if nargin < 6
    option.version = 'FH';
    option.frame = 1;
end
switch option.version
    case 'FH'
%         keyboard
        M = hessianMatrix(H.H,H.s,H.y,H.delta,H.Ginv0,H.i); % preconditioner
        x = x_start;
        r0 = A*(x) - b;
        r = r0;
%         keyboard
        m = M*r0;
        p = m;
%         p = -r0;
        color = dre;
        % ----DEBUG----
        err=[];
        %{%
        if isempty(H.s)
            NOP;
        else       
%             err = norm(r0);
%             q = A*p;
%             alpha = - (p'*q + 1e-30)\(p'*r);
%             x = x_start + alpha * p;
%             r0 = A * x - b;
            p0 = M * r0;
            r = r0;
            p = p0;
        end
        %}
        % ----DEBUG----
    case 'CG'
        M = eye(size(A)); % preconditioner
%         M = hessianMatrix(H.H,H.s,H.y,H.delta,H.Ginv0,H.i); % preconditioner
        x = x_start;
        r0 = A*(x) - b;
        r = r0;
        m = M*r0;
        p = m;
        color = mpg;
        err=[];
    case 'SU'
        M = hessianMatrix(H.H,H.s,H.y,H.delta,H.Ginv0,H.i); % preconditioner
        H_crt = hessianMatrix(H.H);                     % hessian matrix for current frame
        x = x_start;
        r0 = A*(x) - b;
        r = r0;
        m = M*r0;
%         p = m;  % -(H*r);  % #### DEBUG
        p = -r0;
        color = ora;
        % ----DEBUG----
        err=[];
        if isempty(H.s) && isempty(H.l) 
            NOP;
        else            
            err = norm(r0);
            q = A*p;
            alpha = - (p'*q + 1e-30)\(p'*r);
            x = x_start + alpha * p;
            r0 = A * x - b;
            p0 = M * r0;
            r = r0;
            p = p0;
        end
        % ----DEBUG----
    otherwise
        error('check option.version in pncg_Hmfd')
end


A_inv = inv(A); diffz = [];
epsl = 1e-30;
% err= []; 
x_stack = [];residual = r; %zeros(60,1);
% keyboard
for k = 1:numel(b)
    % all intermediat solutions
    x_stack = [x_stack, x];
    %
    err = [err,norm(r)];
    figure(2), plot(err,'Color',color), drawnow, hold on, set(gca,'Yscale','log')
    if k == iter + 1
        break
    end
    if norm(r) < tol
        disp('==> solution found!')
        break
    else
        % ###################
        q = A*p;
        alpha = - (p'*q + epsl)\(p'*r);
        
        s = alpha*p;           % s_i <-- x_i+1 - x_i
        y = alpha*q;           % y_i <-- A*s_i
        switch option.version
            case 'FH'
%                 delta = s - H*y;        % delta_i <-- s_i - H_i*y_i
                delta = s - y;        % delta_i <-- s_i - H_i*y_i
            case 'CG'
                delta = s - y;        % delta_i <-- s_i - H_i*y_i
            case 'SU'
                delta = s - H_crt*y;        % delta_i <-- s_i - H_i*y_i
% % %                 H_crt = plus(H_crt,s,y,delta); % H_i+1 <-- H_i + (update)
                if s'*y > -1e-15
                    H_crt = plus(H_crt,s,y,delta); % H_i+1 <-- H_i + (update)
                else
                    keyboard
                    display 'direction changes too small.'
                    break
                end
            otherwise
                error('check option.version in pncg_Hmfd')
        end
             
        x = x + s;              % x_i+1 <-- x_i - alpha*p_i
        r = r + y;              % r_i+1 <-- r_i - A*alpa*p_i        
       
% % %         H = plus(H,s,y,delta); % H_i+1 <-- H_i + (update)
        if s'*y > -1e-15
            if option.frame > inf
            % eliminate the component in new tuple [s,y,delta] which are
            % already contained in previous spanned Krylov subspace.
                [S,Y,Delta] = purify(H,s,y,delta); 
%                 keyboard
                clear H
                H = hessianMatrix(M.H,S,Y,Delta,(M.i+1));
            else
                H = plus(H,s,y,delta); % H_i+1 <-- H_i + (update)
            end
        else
                    keyboard
            display 'direction changes too small.'
            break
        end
        
        p_1 = p;
        switch option.version
            case 'FH'
                p = H*r;               % p = H_i+1 * r_i+1
            case 'CG'
                p = r + H.*r;
            case 'SU'
                p = H_crt*r;                    
            otherwise
                error('check option.version in pncg_Hmfd')
        end
%         keyboard
        buildH;
        diffz = [diffz, norm(H_mtx - A_inv)/numel(H.H)];
        figure(997); hold on, plot(diffz), drawnow, xlabel('#steps'),ylabel('norm( H - A^{-1} ) / pixel')
% orthogonal      
        residual = [residual,r]; 
        display(sprintf('orth residual: %d', residual(:,end-1)'*residual(:,end)))
% conjugate
        conj = p_1'*A*p
        % ###################
    end
    
end
figure(22), set(gcf,'visible','off'), plot(err,'Color',color), drawnow,hold on, set(gca,'Yscale','log')
% keyboard

clear H_crt
end