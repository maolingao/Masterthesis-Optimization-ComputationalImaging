% plot for blind
startup;

figPath = '/home/gao/Documents/MPI/thesis/article/figure/fig4thesis';

fclk_err    =   figure;
fclk_rerr   =   figure;
fstp_err    =   figure;
fstp_rerr   =   figure;

%% psf
load('dataK_b_pcg.mat');

setnum     =   length(dataK);
for k = 1 : 30
    timeLine    =   dataK{k}.time;
    stepLine    =   1:length(dataK{k}.errs);
    color = dre;
    LineStyle = '-';
    LineWidth = 1;
    
    figure(fclk_err),   
    hData = plot(timeLine,dataK{k}.errs,'Color',color,'lineStyle',LineStyle,'lineWidth',LineWidth); hold on
    figure(fclk_rerr),
    hData = plot(timeLine,dataK{k}.rerrs,'Color',color,'lineStyle',LineStyle,'lineWidth',LineWidth); hold on
    figure(fstp_err),  
    hData = plot(stepLine,dataK{k}.errs,'Color',color,'lineStyle',LineStyle,'lineWidth',LineWidth); hold on
    figure(fstp_rerr),  
    hData = plot(stepLine,dataK{k}.rerrs,'Color',color,'lineStyle',LineStyle,'lineWidth',LineWidth); hold on
end

%%
figtoken = 'b_pcg';

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
%
figure(fstp_err)
hXLabel = xlabel('$\#\text{steps}$', 'Interpreter','Latex');
hYLabel = ylabel('$\|f * x - y\| / pixel$'); 
thisFigure
% xlim([1 4]);
figName = strcat('fstp_err_',figtoken);
figName = strcat(figName,'.tikz');
figname = fullfile(figPath,figName);
printTikz;
%
figure(fstp_rerr)
hXLabel = xlabel('$\#\text{steps}$', 'Interpreter','Latex');
hYLabel = ylabel('$\text{relative error}$'); 
thisFigure
% xlim([1 4]);
figName = strcat('fstp_rerr_',figtoken);
figName = strcat(figName,'.tikz');
figname = fullfile(figPath,figName);
printTikz;
