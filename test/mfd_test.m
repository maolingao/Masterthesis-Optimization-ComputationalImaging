% mfd_test
clear all
%
startup;
localsetup;
figPath = option.figPath;
%
n = 60;
Q = RandomRotation(n); 
u = rand(n,1);
u = clip(u,1,0.1);
u = sort(u,'ascend');
step = 50; u(1:step) = 1e-1*u(1:step)+0.1; u(step+1:end) = u(step+1:end);
% step = 50; u(1:step) = 0.1*ones(step,1); u(step+1:end) = u(step+1:end);
D = diag(u);
A = Q*D*Q';
H_true = Q * diag(1./u) * Q';
% b = randn(n,1);
% x = H_true * b;
tol = 1e-8;
iter = option.iter;

% Identity + H_true
% option.H0fun = @(x) eye(size(x,1))*x;
% option.Wfun = @(x) H_true*x;

% Identity + Identity
% option.H0fun = @(x) eye(size(x,1))*x;
% option.Wfun = @(x) eye(size(x,1))*x;

%% CG & PCG
figure(997), clf
figure(2), clf
figure(22), set(gcf,'visible','off'),clf
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
option.MEMLIM = 60;     % <-- toggle here for different MEMLIM for same problem(ONLY excute this cell!)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%
option.solverMode = 'CG';
option.flag_pa = 0;
step =  [];
R0 = Q(:,step); D0 = diag(1./u(step) - 1);
H_approx = R0 * D0 * R0';
option.H0fun = @(x) x + H_approx*x;

switch option.solverMode
    case 'Greenstadt'
        option.Wfun = option.H0fun;
        if option.flag_pa == 0 || ~isfield(option,'flag_pa')
            option.linestyle = '-';
            H = hessianMatrix(eye(size(A)),[],[],[],[],[],option.Wfun,option.H0fun); % qn
        else
            option.linestyle = '-.';
            H = hessianMatrix(eye(size(A)),[],[],[],R0,D0); % pa
        end
    case 'CG'
        option.Wfun = @(x) H_true*x;
        if option.flag_pa == 0 || ~isfield(option,'flag_pa')
            option.linestyle = '-';
            H = hessianMatrix(eye(size(A)),[],[],[],[],[],option.Wfun,option.H0fun); % qn
        else
            option.linestyle = '-.';
            H = hessianMatrix(eye(size(A)),[],[],[],R0,D0,option.Wfun); % pa
        end
    otherwise
        error('malformed option.solverMode.s')
end

for i = 1: option.numFrame

% for i = 1: option.numFrame
b = rand(n,1);
x_start =   zeros(size(b)); 
% cg solver
% [x_cg,x_cg_seq] = cg(A,b,x_start,tol,iter);
% pcg solver
option.colorIdx = i*10;
[x_pncg,H,residual_pncg,x_pncg_seq]= pncg_Hmfd(A,b,H,x_start,tol,iter,option);
H
% ########### FIGURE #############
% Gram matrix
% figure(101), set(gcf,'visible','off')
% imagesc(log10(abs((H.y'*option.Wfun(H.y))))), colormap gray,  axis image
% set(gca,'xtick',[],'ytick',[]); box on 
% xlabel('column'); ylabel('row'); colorbar('southoutside'); 
% figname = strcat('gm_toy_',option.version,'_', num2str(i), '.tikz');
% figname = fullfile(figPath,figname);
% printTikz;
% ########### MEMLIM #############
%{%
MEMLIM = option.MEMLIM;
lambda = option.MEMSTR;
alpha  = option.EXPOSTR;
% ------------------------------ %
% ### evd
[S,Y,Delta,GInv] = purify(H.s,H.y,H.delta,MEMLIM,lambda);
clear H
H = hessianMatrix(eye(size(A)), S, Y, Delta, [],[],option.Wfun,option.H0fun);
% ------------------------------ %
% ### low rank evd
% option.data.R = H.R;
% option.data.D = H.D;
% option.Wfun = H.Wfun;
% option.H0fun = H.H0fun;
% 
% [R,D] = purify_lowRank(H.s,H.y,H.delta,MEMLIM,H.R,H.D,option);
% 
% H_approx = R * D * R';
% clear H  option.H0fun option.Wfun
% option.H0fun = @(x) x + H_approx*x;
% option.Wfun = @(x) H_true*x;
% H = hessianMatrix(eye(size(A)),[],[],[],R,D,option.Wfun);%,option.H0fun); % qn
%}
end
% H_mtx -> A^{-1}
figure(997)
figname = strcat('H_toy_',option.version,'_', num2str(i), '.tikz');
figname = fullfile(figPath,figname);
printTikz;

% residual
figure(22),set(gcf,'visible','off'); 
% hLegend = legend('classic','probabilistic');
% hLegend = legend('classic','rank10','rank20','rank30','rank40','rank50','rank60');
hLegend = legend('frame 1','frame 2','frame 2','frame 4','frame 5','frame 6','frame 7','frame 8','frame 9','frame 10');
% hTitle = title('H0fun:lra; Wfun:H\_true');
set(hLegend,'Location','northeast')
hYLabel = ylabel('$\|Bx - b\| / pixel$');
hXLabel = xlabel('$\#\text{steps}$');
thisFigure;

figname = sprintf('multiframe_evdG_firstStep_QN_clustered_MEMLIM%d.tikz',option.MEMLIM); figname = fullfile(figPath,figname);
printTikz;
%
%
%
%
% convert eps to pdf
% system(sprintf('/is/ei/mgao/own_scripts/eps2pdf_ule/eps2pdf_ule /is/ei/mgao/Documents/thesis/notes/cgs/fig/*.eps'))
