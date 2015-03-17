% mfd_test
clear all
figure(997), clf
figure(2), clf
figure(22), set(gcf,'visible','off'),clf
%
localsetup;
figPath = option.figPath;
%
n = 60;
% u =-10 * log(rand(n,1));u(1:5) = 100*u(1:5);
u = rand(n,1);  
% u = clip((u*10),10,0); 
% step = 10; u(1:step) = 10*u(1:step); u(step+1:end) = u(step+1:end)./10;
Q = RandomRotation(n); 
D = diag(u);
A = Q*D*Q';
tol = 1e-14;
iter = option.iter;
%% CG & PCG
H = hessianMatrix(eye(size(A)));
for i = 1: option.numFrame

x = rand(n,1);
b = A*x;
x_start =   zeros(size(b)); 
% cg solver
[x_cg] = cg(A,b,x_start,tol,iter);
% pcg solver
[x_pncg,H,residual_pncg]= pncg_Hmfd(A,b,H,x_start,tol,iter,option);
% ########### MEMLIM #############
%{%
MEMLIM = option.MEMLIM;
lambda = option.MEMSTR;
alpha  = option.EXPOSTR;
% ------------------------------ %
% ### evd
% [S,Y,Delta,GInv] = purify(H.s,H.y,H.delta,MEMLIM,lambda);
% clear H
% H = hessianMatrix(eye(size(A)), S, Y, Delta, [],[], GInv);
% ------------------------------ %
% ### low rank evd
keyboard
[R,D] = purify_lowRank(H.s,H.y,H.delta,MEMLIM,H.R,H.D);
clear H
H = hessianMatrix(eye(size(A)), [], [], [], R, D);
%}
% ################################

% ########### FIGURE #############
% Gram matrix
figure(101), set(gcf,'visible','off')
imagesc(log10(abs((H.s'*H.y)))), colormap gray,  axis image off
colorbar('southoutside')
figname = strcat('gm_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(figPath,figname);
print('-depsc2', figname);

end

% H_mtx -> A^{-1}
figure(997)
figname = strcat('H_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(figPath,figname);
print('-depsc2',figname)

% residual
figure(22),set(gcf,'visible','off'); 
hLegend = legend('classic','probabilistic');
hYLabel = ylabel('$\|Bx - b\| / pixel$');
hXLabel = xlabel('$\#steps$');
thisFigure;
figname = 'residualPCGCG.eps'; figname = fullfile(figPath,figname);
print('-depsc2',figname)
%
%
%
%
% convert eps to pdf
% system(sprintf('/is/ei/mgao/own_scripts/eps2pdf_ule/eps2pdf_ule /is/ei/mgao/Documents/thesis/notes/cgs/fig/*.eps'))
