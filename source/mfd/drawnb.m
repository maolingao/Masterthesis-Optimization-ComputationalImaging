% draw residual comparison curves - multiframe non-blind
startup;
figPath = '/home/gao/Documents/MPI/thesis/article/figure/fig4thesis';

fclk_err    =   figure;
fclk_rerr   =   figure;
fstp_err    =   figure;
fstp_rerr   =   figure;
%% 4 legend
load('data_nb_cg.mat');

timeLine    =   data.time;
stepLine    =   1:size(data.errs,1);

figure(fclk_err),   
hData = plot(timeLine(:,1),data.errs(:,1),'Color', mpg, 'linewidth',1); hold on
hData = plot(timeLine(:,1),data.errs(:,1),'Color', dre, 'linewidth',1);
figure(fclk_rerr),
hData = plot(timeLine(:,1),data.rerrs(:,1),'Color', mpg, 'linewidth',1); hold on
hData = plot(timeLine(:,1),data.rerrs(:,1),'Color', dre, 'linewidth',1);
figure(fstp_err),  
hData = plot(stepLine,data.errs(:,1),'Color', mpg, 'linewidth',1); hold on
hData = plot(stepLine,data.errs(:,1),'Color', dre, 'linewidth',1);
figure(fstp_rerr),  
hData = plot(stepLine,data.rerrs(:,1),'Color', mpg, 'linewidth',1); hold on
hData = plot(stepLine,data.rerrs(:,1),'Color', dre, 'linewidth',1);

%% cg
load('data_nb_cg.mat')
timeLine            =   data.time;
for k = 1 : size(data.time,2)-1
    timeLine(:,k+1) = timeLine(:,k+1) + timeLine(end,k);
end
timeLine = vec(timeLine);
stepLine    =   1:numel(data.errs);

figure(fclk_err),   
hData = plot(timeLine,vec(data.errs),'Color', mpg, 'linewidth',1); hold on
figure(fclk_rerr),
hData = plot(timeLine,vec(data.rerrs),'Color', mpg, 'linewidth',1); hold on
figure(fstp_err),  
hData = plot(stepLine,vec(data.errs),'Color', mpg, 'linewidth',1); hold on
figure(fstp_rerr),  
hData = plot(stepLine,vec(data.rerrs),'Color', mpg, 'linewidth',1); hold on

%% pcg
load('data_nb_pcg.mat')
timeLine            =   data.time;
for k = 1 : size(data.time,2)-1
    timeLine(:,k+1) = timeLine(:,k+1) + timeLine(end,k);
end
timeLine = vec(timeLine);
[timeLine,idx,~]           =   unique(timeLine);

stepLine    =   1:length(idx);

figure(fclk_err),   
hData = plot(timeLine,vec(data.errs(idx)),'Color', dre, 'linewidth',1); hold on
figure(fclk_rerr),
hData = plot(timeLine,vec(data.rerrs(idx)),'Color', dre, 'linewidth',1); hold on
figure(fstp_err),  
hData = plot(stepLine,vec(data.errs(idx)),'Color', dre, 'linewidth',1); hold on
figure(fstp_rerr),  
hData = plot(stepLine,vec(data.rerrs(idx)),'Color', dre, 'linewidth',1); hold on
%% 
figure(fclk_err)
hXLabel = xlabel('$\text{time/sec}$', 'Interpreter','Latex');
hYLabel = ylabel('$\|f * x - y\| / pixel$');
hLegend = legend('CG classic', 'CG prob');%, 'R-Lucy', 'Gaussian');
set(hLegend,'location','northeast')
thisFigure
figName = strcat('fclk_err_all','.tikz');
figname = fullfile(figPath,figName);
printTikz;
%
figure(fclk_rerr)
hXLabel = xlabel('$\text{time/sec}$', 'Interpreter','Latex');
hYLabel = ylabel('$\text{relative error}$');
hLegend = legend('CG classic', 'CG prob');%, 'R-Lucy', 'Gaussian');
set(hLegend,'location','northeast')
thisFigure
figName = strcat('fclk_rerr_all','.tikz');
figname = fullfile(figPath,figName);
printTikz;
%
figure(fstp_err)
hXLabel = xlabel('$\#\text{steps}$', 'Interpreter','Latex');
hYLabel = ylabel('$\|f * x - y\| / pixel$');
hLegend = legend('CG classic', 'CG prob');%, 'R-Lucy', 'Gaussian');
set(hLegend,'location','northeast')
thisFigure
figName = strcat('fstp_err_all','.tikz');
figname = fullfile(figPath,figName);
printTikz;
%
figure(fstp_rerr)
hXLabel = xlabel('$\#\text{steps}$', 'Interpreter','Latex');
hYLabel = ylabel('$\text{relative error}$');
hLegend = legend('CG classic', 'CG prob');%, 'R-Lucy', 'Gaussian');
set(hLegend,'location','northeast')
thisFigure
figName = strcat('fstp_rerr_all','.tikz');
figname = fullfile(figPath,figName);
printTikz;