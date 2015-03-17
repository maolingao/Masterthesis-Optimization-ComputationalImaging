function [x,H,residual] = pncg_Hmfd(A,b,H,x_start,tol,iter,option)
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
    iter = 10;
end
if nargin < 7
    option.version = 'FH';
    option.frame = 1;
end

switch option.version
    case 'FH'
        M = hessianMatrix(H.H,H.s,H.y,H.delta,H.R,H.D,H.Ginv0);
        x = x_start;
        r = A*(x) - b;
        p = M*r;
        color = dre;
        
    case 'CG'
        M = eye(size(A));
        x = x_start;
        r = A*(x) - b;
        p = M*r;
        color = mpg;
    otherwise
        error('[pncg_Hmfd.m] : option.version can be either FH or CG.')
end


A_inv = inv(A); 
diffz = []; err=[];
epsl = 1e-30;
residual = r; 

for k = 1:numel(b)
    
    err = [err,norm(r)];
    figure(2), hData = plot(0:length(err)-1,err,'Color',color); drawnow, hold on, set(gca,'Yscale','log')
    thisFigure;
    
    if k == iter + 1
        disp('==> maximal iteration reached!')
        break
    end
    if norm(r) < tol
        disp('==> solution found!')
        break
    else
% ########### MAIN PART #############
        q = A*p;
        alpha = - (p'*q + epsl)\(p'*r);
        
        s = alpha*p;           % s_i <-- x_i+1 - x_i
        y = alpha*q;           % y_i <-- A*s_i
        
        s_length = norm(s,'fro');
        s = s ./ s_length;
        y = y ./ s_length;
        
        switch option.version
            case 'FH'
                delta = s - y;        % delta_i <-- s_i - H_i*y_i
            case 'CG'
                delta = s - y;        % delta_i <-- s_i - H_i*y_i
        end
        
        x = x +  s_length*s;              % x_i+1 <-- x_i - alpha*p_i
        r = r +  s_length*y;              % r_i+1 <-- r_i - A*alpa*p_i        
       
        if abs(s'*y) > 1e-15
            H = plus(H,s,y,delta);        % H_i+1 <-- H_i + (update)
        else
            disp('==> update changes too small!')
            break
        end
        
        p_1 = p;
        switch option.version
            case 'FH'
                p = H*r;               % p = H_i+1 * r_i+1
            case 'CG'
                p = r + H.*r;
        end
        
% ########### DIAGNOSTIC #############
% build H_mtx explicitly --> A^{-1}
        buildH;
        diffz = [diffz, norm(H_mtx - A_inv) / numel(H.H)];
        figure(997); 
        hData = plot(diffz,'b'); drawnow, hold on;
        hXLabel = xlabel('$\#steps$'); hYLabel = ylabel('$\| H - B^{-1} \| / pixel$');
        thisFigure;
% orthogonal      
        residual = [residual,r]; 
        display(sprintf('orth residual: %d', residual(:,end-1)'*residual(:,end)))
% conjugate
        conj = p_1'*A*p
    end
    
end
% residual curve for latex
figure(22), set(gcf,'visible','off'), 
hData = plot(0:length(err)-1,err,'Color',color); 
drawnow, hold on;
set(gca,'Yscale','log');
thisFigure;

end