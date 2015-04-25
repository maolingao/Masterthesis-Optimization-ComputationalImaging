% draw residual comparison curves - multiframe non-blind
startup;
figPath = '/home/gao/Documents/MPI/thesis/article/figure/fig4thesis';

fclk_err    =   figure;
fclk_rerr   =   figure;
fstp_err    =   figure;
fstp_rerr   =   figure;
% figure,
% fclk_err    =   figure(221);
% fclk_rerr   =   figure(222);
% fstp_err    =   figure(223);
% fstp_rerr   =   figure(224);
%% 4 legend
load('data_nb_cg.mat');

timeLine    =   data.time;
stepLine    =   1:size(data.errs,1);

figure(fclk_err),   
hData = loglog(timeLine(:,1),data.errs(:,1),'Color', mpg, 'linewidth',1); hold on
hData = loglog(timeLine(:,1),data.errs(:,1),'Color', dre, 'linewidth',1); hold on
hData = loglog(timeLine(:,1),data.errs(:,1),'Color', ora, 'linewidth',1); hold on
hData = loglog(timeLine(:,1),data.errs(:,1),'Color', blu, 'linewidth',1); 
figure(fclk_rerr),
hData = loglog(timeLine(:,1),data.rerrs(:,1),'Color', mpg, 'linewidth',1); hold on
hData = loglog(timeLine(:,1),data.rerrs(:,1),'Color', dre, 'linewidth',1); hold on
hData = loglog(timeLine(:,1),data.rerrs(:,1),'Color', ora, 'linewidth',1); hold on
hData = loglog(timeLine(:,1),data.rerrs(:,1),'Color', blu, 'linewidth',1);
figure(fstp_err),  
hData = loglog(stepLine,data.errs(:,1),'Color', mpg, 'linewidth',1); hold on
hData = loglog(stepLine,data.errs(:,1),'Color', dre, 'linewidth',1); hold on
hData = loglog(stepLine,data.errs(:,1),'Color', ora, 'linewidth',1); hold on
hData = loglog(stepLine,data.errs(:,1),'Color', blu, 'linewidth',1);
figure(fstp_rerr),  
hData = loglog(stepLine,data.rerrs(:,1),'Color', mpg, 'linewidth',1); hold on
hData = loglog(stepLine,data.rerrs(:,1),'Color', dre, 'linewidth',1); hold on
hData = loglog(stepLine,data.rerrs(:,1),'Color', ora, 'linewidth',1); hold on
hData = loglog(stepLine,data.rerrs(:,1),'Color', blu, 'linewidth',1);

%% cg
load('data_nb_cg.mat');

timeLine    =   data.time;
stepLine    =   1:size(data.errs,1);

    figure(fclk_err),   hData = loglog(timeLine,data.errs,'Color', mpg, 'linewidth',1); hold on, thisFigure
    figure(fclk_rerr),  hData = loglog(timeLine,data.rerrs,'Color', mpg, 'linewidth',1); hold on, thisFigure
    figure(fstp_err),   hData = loglog(stepLine,data.errs,'Color', mpg, 'linewidth',1); hold on, thisFigure
    figure(fstp_rerr),  hData = loglog(stepLine,data.rerrs,'Color', mpg, 'linewidth',1); hold on, thisFigure

%% pcg
clear data
load('data_nb_pcg.mat');

timeLine    =   data.time;
stepLine    =   1:size(data.errs,1);

    figure(fclk_err),   hData = loglog(timeLine,data.errs,'Color', dre, 'linewidth',1); hold on, thisFigure
    figure(fclk_rerr),  hData = loglog(timeLine,data.rerrs,'Color', dre, 'linewidth',1); hold on, thisFigure
    figure(fstp_err),   hData = loglog(stepLine,data.errs,'Color', dre, 'linewidth',1); hold on, thisFigure
    figure(fstp_rerr),  hData = loglog(stepLine,data.rerrs,'Color', dre, 'linewidth',1); hold on, thisFigure

% %% rl
% clear data
% load('data_nb_rl.mat');
% 
% timeLine    =   data.time;
% stepLine    =   1:size(data.errs,1);
% 
% figure(fclk_err),   hData = loglog(timeLine,data.errs,'Color', ora, 'linewidth',1); hold on, thisFigure
% figure(fclk_rerr),  hData = loglog(timeLine,data.rerrs,'Color', ora, 'linewidth',1); hold on, thisFigure
% figure(fstp_err),   hData = loglog(stepLine,data.errs,'Color', ora, 'linewidth',1); hold on, thisFigure
% figure(fstp_rerr),  hData = loglog(stepLine,data.rerrs,'Color', ora, 'linewidth',1); hold on, thisFigure
% 
% %% gaussian
% clear data
% load('data_nb_gaussian.mat');
% 
% timeLine    =   data.time;
% stepLine    =   1:size(data.errs,1);
% 
% figure(fclk_err),   hData = loglog(timeLine,data.errs,'Color', blu, 'linewidth',1); hold on, thisFigure
% figure(fclk_rerr),  hData = loglog(timeLine,data.rerrs,'Color', blu, 'linewidth',1); hold on, thisFigure
% figure(fstp_err),   hData = loglog(stepLine,data.errs,'Color', blu, 'linewidth',1); hold on, thisFigure
% figure(fstp_rerr),  hData = loglog(stepLine,data.rerrs,'Color', blu, 'linewidth',1); hold on, thisFigure


%%
figure(fclk_err)
hXLabel = xlabel('$\text{time/sec}$', 'Interpreter','Latex');
hYLabel = ylabel('$\|f * x - y\| / pixel$');
hLegend = legend('CG classic', 'CG prob');%, 'R-Lucy', 'Gaussian');
set(hLegend,'location','southwest')
thisFigure
figName = strcat('fclk_err_gt','.tikz');
figname = fullfile(figPath,figName);
printTikz;
%
figure(fclk_rerr)
hXLabel = xlabel('$\text{time/sec}$', 'Interpreter','Latex');
hYLabel = ylabel('$\text{relative error}$');
hLegend = legend('CG classic', 'CG prob');%, 'R-Lucy', 'Gaussian');
set(hLegend,'location','southwest')
thisFigure
figName = strcat('fclk_rerr_gt','.tikz');
figname = fullfile(figPath,figName);
printTikz;
%
figure(fstp_err)
hXLabel = xlabel('$\#\text{steps}$', 'Interpreter','Latex');
hYLabel = ylabel('$\|f * x - y\| / pixel$');
hLegend = legend('CG classic', 'CG prob');%, 'R-Lucy', 'Gaussian');
set(hLegend,'location','southwest')
thisFigure
figName = strcat('fstp_err_gt','.tikz');
figname = fullfile(figPath,figName);
printTikz;
%
figure(fstp_rerr)
hXLabel = xlabel('$\#\text{steps}$', 'Interpreter','Latex');
hYLabel = ylabel('$\text{relative error}$');
hLegend = legend('CG classic', 'CG prob');%, 'R-Lucy', 'Gaussian');
set(hLegend,'location','southwest')
thisFigure
figName = strcat('fstp_rerr_gt','.tikz');
figname = fullfile(figPath,figName);
printTikz;