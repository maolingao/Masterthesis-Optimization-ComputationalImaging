% plot for blind % plot for blind

startup;

figPath = '/home/gao/Documents/MPI/thesis/article/figure/fig4thesis';

fclk_err    =   figure;
fclk_rerr   =   figure;
fstp_err    =   figure;
fstp_rerr   =   figure;
%% estimated g.t. all frame error - pcg
load('dataN_b_pcg.mat');

load('dataK_b_pcg.mat');
timeLine    =   calcTimeLine(dataK);

stepLine    =   1:length(dataN.errs);

figure(fclk_err),   
hData = plot(timeLine,dataN.errs,'Color', dre, 'linewidth',1); hold on
figure(fclk_rerr),
hData = plot(timeLine,dataN.rerrs,'Color', dre, 'linewidth',1); hold on
figure(fstp_err),  
hData = plot(stepLine,dataN.errs,'Color', dre, 'linewidth',1); hold on
figure(fstp_rerr),  
hData = plot(stepLine,dataN.rerrs,'Color', dre, 'linewidth',1); hold on
%% estimated g.t. all frame error - cg
load('dataN_b_cg.mat');

load('dataK_b_cg.mat');
timeLine    =   calcTimeLine(dataK);

stepLine    =   1:length(dataN.errs);

figure(fclk_err),   
hData = plot(timeLine,dataN.errs,'Color', mpg, 'linewidth',1); hold on
figure(fclk_rerr),
hData = plot(timeLine,dataN.rerrs,'Color', mpg, 'linewidth',1); hold on
figure(fstp_err),  
hData = plot(stepLine,dataN.errs,'Color', mpg, 'linewidth',1); hold on
figure(fstp_rerr),  
hData = plot(stepLine,dataN.rerrs,'Color', mpg, 'linewidth',1); hold on
%%
figtoken = 'b_gt_cg_pcg';

figure(fclk_err)
hXLabel = xlabel('$\text{time/sec}$', 'Interpreter','Latex');
hYLabel = ylabel('$\|f * x - y\| / pixel$'); 
thisFigure
figName = strcat('fclk_err_',figtoken);
figName = strcat(figName,'.tikz');
figname = fullfile(figPath,figName);
printTikz;
%
figure(fclk_rerr)
hXLabel = xlabel('$\text{time/sec}$', 'Interpreter','Latex');
hYLabel = ylabel('$\text{relative error}$'); 
thisFigure
figName = strcat('fclk_rerr_',figtoken);
figName = strcat(figName,'.tikz');
figname = fullfile(figPath,figName);
printTikz;

figure(fstp_err)
hXLabel = xlabel('$\#\text{frames}$', 'Interpreter','Latex');
hYLabel = ylabel('$\|f * x - y\| / pixel$'); 
thisFigure
figName = strcat('fstp_err_',figtoken);
figName = strcat(figName,'.tikz');
figname = fullfile(figPath,figName);
printTikz;
%
figure(fstp_rerr)
hXLabel = xlabel('$\#\text{frames}$', 'Interpreter','Latex');
hYLabel = ylabel('$\text{relative error}$'); 
thisFigure
figName = strcat('fstp_rerr_',figtoken);
figName = strcat(figName,'.tikz');
figname = fullfile(figPath,figName);
printTikz;