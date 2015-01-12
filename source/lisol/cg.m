function [x] = cg(A,b,x_start,tol,iter)

% example
%{
  n = 6000;
  m = 8000;
  A = randn(n,m);
  A = A * A';
  b = randn(n,1);
  tic, x = cg(A,b); toc
  norm(A*x-b)
%}
startup;
if nargin < 3
    x_start = b;
end

if nargin < 4
    tol = 10^-10;
end
if nargin < 5
    iter = 100;
end
x = x_start;
r = A*x - b;
p = -r;
epsl = 1e-30; % numerical stability

err= [];

for k = 1:numel(b)
    err = [err,norm(r)];
    figure(2), plot(err,'Color',blu), drawnow,hold on, set(gca,'Yscale','log');
    if k == iter + 1
        break
    end
    if norm(r) < tol
        disp('==> solution found!')
        break
    else
        q = A*p; % A*p register
        alpha = ((p'*q) + epsl)\(r'*r);
        x = x + alpha*p;
        r_1 = r;
        r = r_1 + alpha*q;
        beta = ((r_1'*r_1) + epsl)\(r'*r);
        p_1 = p;
        p = -r + beta*p_1;
        conj = p_1'*A*p
    end
    
end

figure(22), set(gcf,'visible','off'), plot(err,'Color',blu), drawnow,hold on, set(gca,'Yscale','log')


end