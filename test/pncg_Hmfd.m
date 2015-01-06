function [x,H,x_stack] = pncg_Hmfd(A,b,H,x_start,tol,iter,option)
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
end
switch option.version
    case 'FH'
        M = eye(size(A)); % preconditioner
        x = x_start;
        r0 = b -A*(x);
        r = r0;
        m = M*r0;
        p = m;
        color = dre;
    case 'CG'
        M = eye(size(A)); % preconditioner
        x = x_start;
        r0 = b -A*(x);
        r = r0;
        m = M*r0;
        p = m;
        color = mpg;
    case 'SU'
        M = hessianMatrix(H.H,H.l,H.r,H.s,H.y,H.delta); % preconditioner
        H_crt = hessianMatrix(H.H);                     % hessian matrix for current frame
        x = x_start;
        r0 = b - A*x;
        r = r0;
        m = M*r0;
        p = m ; % -(H*r);  % #### DEBUG
        color = ora;
    otherwise
        error('check option.version in pncg_Hmfd')
end


A_inv = inv(A); diffz = [];
epsl = 1e-30;
err= []; x_stack = [];
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
                delta = s - H*y;        % delta_i <-- s_i - H_i*y_i
            case 'CG'
                delta = s - y;        % delta_i <-- s_i - H_i*y_i
            case 'SU'
                delta = s - H_crt*y;        % delta_i <-- s_i - H_i*y_i
                H_crt = plus(H_crt,s,y,delta); % H_i+1 <-- H_i + (update)
            otherwise
                error('check option.version in pncg_Hmfd')
        end
             
        x = x + s;              % x_i+1 <-- x_i - alpha*p_i
        r = r + y;              % r_i+1 <-- r_i - A*alpa*p_i        
%         m = M * r;              % new preconditioned system, define m = M * r
        
%         den = s'*y;
%         H = H + (den + epsl)\(delta*s') + (den + epsl)\(s*delta') - ((den)^2 + epsl)\(s*delta'*y*s');% H_i+1 <-- H_i + (update)
        H = plus(H,s,y,delta); % H_i+1 <-- H_i + (update)
        
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

%         p = r + (den + epsl)\(delta*(s'*r)) + (den + epsl)\(s*(delta'*r)) - ((den)^2 + epsl)\(s*delta'*y*(s'*r));
%         p = r + H.*r;
        
        buildH;
        diffz = [diffz, norm(H_mtx - A_inv)/numel(H.H)];
%         diffz = [diffz, norm(H - A_inv)/numel(H)];
        figure(997); hold on, plot(diffz), drawnow, xlabel('#steps'),ylabel('norm( H - A^{-1} ) / pixel')
% orthogonal
        orth = p_1'*A*p
        % ###################
    end
    
end
figure(22), set(gcf,'visible','off'), plot(err,'Color',color), drawnow,hold on, set(gca,'Yscale','log')
keyboard

end