% draw regularization 

startup;
figPath = '/home/gao/Documents/MPI/thesis/article/figure/fig4thesis';

%%
option.mode = 'eta';
switch option.mode
    case 'SNR'
      load('data_snrDiff_eta0.mat');
      figtoken = 'snrDiff_eta0';
    case 'eta'
      load('data_snrFix_etaDiff.mat');
      figtoken = 'snrFix_etaDiff';
end

figure
fclk_err    =   figure;
fclk_rerr   =   figure;
fstp_err    =   figure;
fstp_rerr   =   figure;

setnum = length(data_snr);

for k = 1 : setnum
    timeLine    =   data_snr{k}.time;
    stepLine    =   1:length(data_snr{k}.errs);
    
    switch k
        case 1
            color = dre; % set color and linestyle of hData
            LineStyle = '-'; % set color and linestyle of hData
            LineWidth = 2;
        case 2
            color = ora; % set color and linestyle of hData
            LineStyle = ':'; % set color and linestyle of hData
            LineWidth = 2;
        case 3
            color = blu; % set color and linestyle of hData
            LineStyle = '-.'; % set color and linestyle of hData
            LineWidth = 2;
        case 4
            color = gra; % set color and linestyle of hData
            LineStyle = '--'; % set color and linestyle of hData
            LineWidth = 2;
        case 5
            color = mpg; % set color and linestyle of hData
            LineStyle = '-.'; % set color and linestyle of hData
            LineWidth = 2;
        case 6
            color = bla; % set color and linestyle of hData
            LineStyle = '-'; % set color and linestyle of hData
            LineWidth = 2;
        otherwise
            error('[drawTikv.m] : malformed eta.')
    end
    
    figure(fclk_err),   
    hData = plot(timeLine,data_snr{k}.errs,'Color',color,'lineStyle',LineStyle,'lineWidth',LineWidth); hold on
    figure(fclk_rerr),
    hData = plot(timeLine,data_snr{k}.rerrs,'Color',color,'lineStyle',LineStyle,'lineWidth',LineWidth); hold on
    figure(fstp_err),  
    hData = plot(stepLine,data_snr{k}.errs,'Color',color,'lineStyle',LineStyle,'lineWidth',LineWidth); hold on
    figure(fstp_rerr),  
    hData = plot(stepLine,data_snr{k}.rerrs,'Color',color,'lineStyle',LineStyle,'lineWidth',LineWidth); hold on
    
end

%% 
    etaVec = [0, 1e-4, 1e-3, 1e-2, 1e-1, 1e0];
    SNRVec = 10:10:60;
    switch option.mode
        case 'SNR'
            trgVec = SNRVec;
            for k = 1 : length(trgVec)
                content{k} = strcat('SNR', sprintf(' %gdB',trgVec(k)));
            end
        case 'eta'
            trgVec = etaVec;
            for k = 1 : length(trgVec)
                content{k} = strcat('$\eta$', sprintf(' %g',trgVec(k)));
            end
        otherwise
            content = [];
    end

figure(fclk_err)
hXLabel = xlabel('$\text{time/sec}$', 'Interpreter','Latex');
hYLabel = ylabel('$\|f * x - y\| / pixel$'); 
hLegend = legend(content); 
set(hLegend,'Interpreter','Latex','Location','northeast');
thisFigure
figName = strcat('fclk_err_',figtoken);
figName = strcat(figName,'.tikz');
figname = fullfile(figPath,figName);
printTikz;

%
figure(fclk_rerr)
hXLabel = xlabel('$\text{time/sec}$', 'Interpreter','Latex');
hYLabel = ylabel('$\text{relative error}$'); 
hLegend = legend(content); 
set(hLegend,'Interpreter','Latex','Location','northeast');
thisFigure
figName = strcat('fclk_rerr_',figtoken);
figName = strcat(figName,'.tikz');
figname = fullfile(figPath,figName);
printTikz;
ylim([0 1]);
figName = strcat('fclk_rerr_',figtoken);
figName = strcat(figName,'_zoom');
figName = strcat(figName,'.tikz');
figname = fullfile(figPath,figName);
printTikz;
%
figure(fstp_err)
hXLabel = xlabel('$\#\text{steps}$', 'Interpreter','Latex');
hYLabel = ylabel('$\|f * x - y\| / pixel$'); 
hLegend = legend(content); 
set(hLegend,'Interpreter','Latex','Location','northeast');
thisFigure
figName = strcat('fstp_err_',figtoken);
figName = strcat(figName,'.tikz');
figname = fullfile(figPath,figName);
printTikz;
%
figure(fstp_rerr)
hXLabel = xlabel('$\#\text{steps}$', 'Interpreter','Latex');
hYLabel = ylabel('$\text{relative error}$'); 
hLegend = legend(content); 
set(hLegend,'Interpreter','Latex','Location','northeast');
thisFigure
figName = strcat('fstp_rerr_',figtoken);
figName = strcat(figName,'.tikz');
figname = fullfile(figPath,figName);
printTikz;
ylim([0 1]);
figName = strcat('fstp_rerr_',figtoken);
figName = strcat(figName,'_zoom');
figName = strcat(figName,'.tikz');
figname = fullfile(figPath,figName);
printTikz;