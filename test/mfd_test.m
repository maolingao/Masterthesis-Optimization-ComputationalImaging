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
u = rand(n,1);  u = clip((u*10),10,0); 
step = 10; u(1:step) = 10*u(1:step); u(step+1:end) = u(step+1:end)./10;
Q = RandomRotation(n); 
D = diag(u);
A = Q*D*Q';
A_ori= A;
tol = 1e-14;
iter = option.iter;
%% precondition
H = hessianMatrix(eye(size(A)));
for i = 1: option.numFrame

x = rand(n,1);
b = A_ori*x;
% b = b./ max(u);
x_start =   zeros(size(b)); %b; %

[x_cg] = cg(A,b,x_start,tol,iter);
%
if i == 1
    option.flag = 0;
else
    option.flag = 1;
end
option.frame = i;
[x_pncg,H,~,residual_pncg]= pncg_Hmfd(A,b,H,x_start,tol,iter,option);
% ########### MEMLIM #############
%{%
MEMLIM = option.MEMLIM ;% size(H.s,2);
lambda = option.memoryStrength;
% keyboard
[S,Y,Delta,GInv] = purify(H.s,H.y,H.delta,H.Ginv0,MEMLIM,lambda);
% keyboard
clear H
H = hessianMatrix(eye(size(A)), S, Y, Delta, [],[], GInv);
%}
% ################################

% Gram matrix
figure(101), imagesc(log10(abs((H.s'*H.y)))), colormap gray,  axis image off
colorbar('southoutside')
figname = strcat('gm_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(figPath,figname);
print('-depsc2', figname);

% #########
% A*H_mtx -> asymptotic identity matrix
buildH;
figure(1), imagesc(log10(abs(A*H_mtx))), colormap gray, axis image off 
colorbar('southoutside')
figname = strcat('HA_toy_',option.version,'_', num2str(i), 'frames.eps');
figname = fullfile(figPath,figname);
print('-depsc2',figname);
end

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
