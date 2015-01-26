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
u = rand(n,1); % u(1:10) = 10*u(1:10);
Q = RandomRotation(n); 
D = diag(u);
A = Q*D*Q';
% figure(5), imagesc(A),colormap
% H  = eye(size(A));
H_CG = hessianMatrix(eye(size(A)));
H_FH = hessianMatrix(eye(size(A)));

%
tol = 1e-14;
iter = 30;
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
    figname = fullfile(figPath,figname);
    print('-depsc2',figname)
    % Gram matrix
    figure(101), imagesc(log10(abs((H_uni.s'*H_uni.y)))), colormap gray,  axis equal 
    colorbar('southoutside')
    figname = strcat('gm_toy_',option.version,'_', num2str(i), '.eps');
    figname = fullfile(figPath,figname);
    print('-depsc2',figname)
    %
    clear H_uni
end

% residual
figure(2)
figname = 'residual.eps';
figname = fullfile(figPath,figname);
print('-depsc2', figname);
% H
figure(997)
figname = strcat('H_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(figPath,figname);
print('-depsc2',figname)
keyboard
%}
%% precondition
%{%
H_SU = hessianMatrix(eye(size(A)));
residual_pncg_allframe = [];
for i = 1: 10

x = rand(n,1);
b = A*x;

if i == 1
    x_start =   zeros(size(b)); %b; %
    r_origin = A * x_start - b;
    b_orth = b;
else
    [b_orth, alpha] = split(b,b_1);
    x_start =   zeros(size(b)); % H_SU * (b_orth + r_origin) ; % % # star #
%     keyboard
end

[x_cg] = cg(A,b_orth,x_start,tol,iter);
if i == 1
    NOP;
else
    [x_cg] = assemble(x_cg, x_cg_1, alpha);
end
x_cg_1 = x_cg;

%
option.version = 'FH';
if i == 1
    option.flag = 0;
else
    option.flag = 1;
end
option.frame = i;
% keyboard
[x_pncg,H_SU,~,residual_pncg]= pncg_Hmfd(A,b_orth,H_SU,x_start,tol,iter,option);
residual_pncg_allframe = [residual_pncg_allframe,residual_pncg];
if i == 1
    NOP;
else
    [x_pncg] = assemble(x_pncg, x_pncg_1, alpha);
end
x_pncg_1 = x_pncg;
b_1 = b;
% ################################
%{%
MEMLIM = 30;% size(H_SU.s,2);
% keyboard
[S,Y,Delta,GInv] = purify(H_SU.s,H_SU.y,H_SU.delta,H_SU.Ginv0,MEMLIM);
% keyboard
clear H_SU
H_SU = hessianMatrix(eye(size(A)), S, Y, Delta, GInv, size(S,2)+1);
max(max(S'*Y - diag(1./diag(H_SU.Ginv0))))
% max(max(diag(1./diag(S'*Y)) - (H_SU.Ginv0))))
%}
% ################################
% keyboard
% H_SU_up = updateH(H_SU,10);

% Gram matrix
figure(101), imagesc(log10(abs((H_SU.s'*H_SU.y)))), colormap gray,  axis image 
colorbar('southoutside')
figname = strcat('gm_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(figPath,figname);
print('-depsc2', figname);
% keyboard
% % % residual matrix
% % figure(102), imagesc(log10(abs((residual_pncg_allframe'*H_SU.s)))), colormap gray,  axis equal 
% % colorbar('southoutside')
% % figname = strcat('rp_toy_',option.version,'_', num2str(i), '.eps');
% % figname = fullfile(figPath,figname);
% % print('-depsc2', figname);

% keyboard
% clear H_SU
% H_SU = H_SU_up;
% clear H_SU_up

% #########
% residual
figure(22),set(gcf,'visible','off'),legend('Nocedal','pncg-SU')
ylabel('$\|Fx - y\| / pixel$','Interpreter','Latex')
xlabel('$\#steps$','Interpreter','Latex')
figname = 'residualSUNocedal.eps';
figname = fullfile(figPath,figname);
print('-depsc2', figname);
% A*H_mtx -> asymptotic identity matrix
clear H, H = H_SU;
buildH;
figure(1), imagesc(log10(abs(A*H_mtx))), colormap gray, axis equal 
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
figure(22),set(gcf,'visible','off'),legend('Nocedal','pncg')
figname = 'residualCGFHNocedal.eps';
figname = fullfile(figPath,figname);
print('-depsc2',figname)
keyboard
%}
%% single and multi frame test
for i = 1: 3

x = rand(n,1);
b = A*x;
x_start = b; % zeros(size(b));  % only allow zero as initial guess
% [x_sol,H]= pncg_Hmfd(A,b,H);
x_cg = cg(A,b,x_start,tol,iter);
%
option.flag = 0;
option.version = 'CG';
[x_pncg,H_CG]= pncg_Hmfd(A,b,H_CG,x_start,tol,iter,option);
% #########
% residual
figure(22),set(gcf,'visible','off'),legend('Nocedal','pncg-exact')
figname = 'residualCGNocedal.eps';
figname = fullfile(figPath,figname);
print('-depsc2', figname);
% A*H_mtx -> asymptotic identity matrix
clear H, H = H_CG;
buildH;
figure(1), imagesc(log10(abs(A*H_mtx))), colormap gray, axis equal 
colorbar('southoutside')
figname = strcat('HA_toy_',option.version,'_', num2str(i), 'frames.eps');
figname = fullfile(figPath,figname);
print('-depsc2',figname)
% Gram matrix
figure(101), imagesc(log10(abs((H_CG.s'*H_CG.y)))), colormap gray,  axis equal 
colorbar('southoutside')
figname = strcat('gm_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(figPath,figname);
print('-depsc2',figname)
%
option.version = 'FH';
[x_pncg,H_FH]= pncg_Hmfd(A,b,H_FH,x_start,tol,iter,option);
% A*H_mtx -> asymptotic identity matrix
clear H, H = H_FH;
buildH;
figure(11), imagesc(log10(abs(A*H_mtx))), colormap gray, axis equal 
colorbar('southoutside')
figname = strcat('HA_toy_',option.version,'_', num2str(i), 'frames.eps');
figname = fullfile(figPath,figname);
print('-depsc2',figname)
% Gram matrix
figure(1001), imagesc(log10(abs((H_FH.s'*H_FH.y)))), colormap gray,  axis equal 
colorbar('southoutside')
figname = strcat('gm_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(figPath,figname);
print('-depsc2',figname)

end
figure(997)
figname = strcat('H_toy_',option.version,'_', num2str(i), '.eps');
figname = fullfile(figPath,figname);
print('-depsc2',figname)
% residual
figure(22),set(gcf,'visible','off'),legend('Nocedal','pncg-exact','pncg-FH')
figname = 'residualCGFHNocedal.eps';
figname = fullfile(figPath,figname);
print('-depsc2',figname)

% convert eps to pdf
% system(sprintf('/is/ei/mgao/own_scripts/eps2pdf_ule/eps2pdf_ule /is/ei/mgao/Documents/thesis/notes/cgs/fig/*.eps'))
