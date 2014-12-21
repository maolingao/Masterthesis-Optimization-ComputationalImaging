% mfd_test
clear all
figure(997), clf
figure(2), clf
figure(22), set(gcf,'visible','off'),clf
pathName = '/is/ei/mgao/Documents/thesis/notes/cgs/fig/';
%
n = 60;
% u =-10 * log(rand(n,1));u(1:5) = 100*u(1:5);
u = rand(n,1);
Q = RandomRotation(n); 
D = diag(u);
A = Q*D*Q';
% figure(5), imagesc(A),colormap
% H  = eye(size(A));
H_CG = hessianMatrix(eye(size(A)));
H_FH = hessianMatrix(eye(size(A)));

%
tol = 1e-14;
iter = 60;
% option.version = 'CG';

%% statistical test
%{
for i = 1: 10
    %    

%     u =-10 * log(rand(n,1));u(1:5) = 100*u(1:5);
    u = 1*rand(n,1); u(u>0.3) = 1;
    Q = RandomRotation(n); 
    D = diag(u);
    A = Q*D*Q';
    %
    H_uni = hessianMatrix(eye(size(A)));
    x = rand(n,1);
    b = A*x;
    x_start = b;
    
%     x_cg = cg(A,b,tol,iter);
    option.version = 'CG';
    [x_sol,H_uni]= pncg_Hmfd(A,b,H_uni,x_start,tol,iter,option);
    
    % A*H_mtx -> asymptotic identity matrix
    clear H, H = H_uni;
    buildH;
    figure(1), imagesc(log10(abs(A*H_mtx))), colormap gray, axis equal 
    colorbar('southoutside')
    figname = strcat('HA_toy_',option.version,'_', num2str(i), 'frames.eps');
    figname = fullfile(pathName,figname);
    print('-depsc2',figname)
    % Gram matrix
    figure(101), imagesc(log10(abs((H_uni.s'*H_uni.y)))), colormap gray,  axis equal 
    colorbar('southoutside')
    figname = strcat('gm_toy_',option.version,'_', num2str(i), '.eps');
    figname = fullfile(pathName,figname);
    print('-depsc2',figname)
    %
    clear H_uni
end

% residual
figure(2)
figname = 'residual.eps';
figname = fullfile(pathName,figname);
print('-depsc2', figname);
% H
figure(997)
figname = strcat('H_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(pathName,figname);
print('-depsc2',figname)
keyboard
%}
%% precondition
%{%
H_SU = hessianMatrix(eye(size(A)));
for i = 1: 5

x = rand(n,1);
b = A*x;
x_start = b; %zeros(size(b));  
% [x_sol,H]= pncg_Hmfd(A,b,H);
x_cg = cg(A,b,tol,iter);
%
option.version = 'SU';
if i == 1
    option.flag = 0;
else
    option.flag = 1;
end
[x_sol,H_SU]= pncg_Hmfd(A,b,H_SU,x_start,tol,iter,option);
% #########
% residual
figure(22),set(gcf,'visible','off'),legend('Nocedal','pncg-SU')
figname = 'residualSUNocedal.eps';
figname = fullfile(pathName,figname);
print('-depsc2', figname);
% A*H_mtx -> asymptotic identity matrix
clear H, H = H_SU;
buildH;
figure(1), imagesc(log10(abs(A*H_mtx))), colormap gray, axis equal 
colorbar('southoutside')
figname = strcat('HA_toy_',option.version,'_', num2str(i), 'frames.eps');
figname = fullfile(pathName,figname);
print('-depsc2',figname)
% Gram matrix
figure(101), imagesc(log10(abs((H_SU.s'*H_SU.y)))), colormap gray,  axis equal 
colorbar('southoutside')
figname = strcat('gm_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(pathName,figname);
print('-depsc2',figname)
end

figure(997)
figname = strcat('H_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(pathName,figname);
print('-depsc2',figname)
% residual
figure(22),set(gcf,'visible','off'),legend('Nocedal','pncg-SU')
figname = 'residualCGFHNocedal.eps';
figname = fullfile(pathName,figname);
print('-depsc2',figname)
keyboard
%}
%% single and multi frame test
for i = 1: 3

x = rand(n,1);
b = A*x;
x_start = b; % zeros(size(b));  % only allow zero as initial guess
% [x_sol,H]= pncg_Hmfd(A,b,H);
x_cg = cg(A,b,tol,iter);
%
option.flag = 0;
option.version = 'CG';
[x_sol,H_CG]= pncg_Hmfd(A,b,H_CG,x_start,tol,iter,option);
% #########
% residual
figure(22),set(gcf,'visible','off'),legend('Nocedal','pncg-exact')
figname = 'residualCGNocedal.eps';
figname = fullfile(pathName,figname);
print('-depsc2', figname);
% A*H_mtx -> asymptotic identity matrix
clear H, H = H_CG;
buildH;
figure(1), imagesc(log10(abs(A*H_mtx))), colormap gray, axis equal 
colorbar('southoutside')
figname = strcat('HA_toy_',option.version,'_', num2str(i), 'frames.eps');
figname = fullfile(pathName,figname);
print('-depsc2',figname)
% Gram matrix
figure(101), imagesc(log10(abs((H_CG.s'*H_CG.y)))), colormap gray,  axis equal 
colorbar('southoutside')
figname = strcat('gm_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(pathName,figname);
print('-depsc2',figname)
%
option.version = 'FH';
[x_sol,H_FH]= pncg_Hmfd(A,b,H_FH,x_start,tol,iter,option);
% A*H_mtx -> asymptotic identity matrix
clear H, H = H_FH;
buildH;
figure(11), imagesc(log10(abs(A*H_mtx))), colormap gray, axis equal 
colorbar('southoutside')
figname = strcat('HA_toy_',option.version,'_', num2str(i), 'frames.eps');
figname = fullfile(pathName,figname);
print('-depsc2',figname)
% Gram matrix
figure(1001), imagesc(log10(abs((H_FH.s'*H_FH.y)))), colormap gray,  axis equal 
colorbar('southoutside')
figname = strcat('gm_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(pathName,figname);
print('-depsc2',figname)

end
figure(997)
figname = strcat('H_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(pathName,figname);
print('-depsc2',figname)
% residual
figure(22),set(gcf,'visible','off'),legend('Nocedal','pncg-exact','pncg-FH')
figname = 'residualCGFHNocedal.eps';
figname = fullfile(pathName,figname);
print('-depsc2',figname)

% convert eps to pdf
% system(sprintf('/is/ei/mgao/own_scripts/eps2pdf_ule/eps2pdf_ule /is/ei/mgao/Documents/thesis/notes/cgs/fig/*.eps'))
