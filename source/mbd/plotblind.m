%% s_snr5_etaf1e-4_etaxDiff
startup;

load('s_snr5_etaf1e-4_etaxDiff');

figure(1)
plot(etax0,'color',dre,'lineWidth',2,'lineStyle','-'), hold on
plot(etax1en1,'color',ora,'lineWidth',2,'lineStyle',':'), hold on
plot(etax1en2,'color',blu,'lineWidth',2,'lineStyle','-.'), hold on
plot(etax1en3,'color',gra,'lineWidth',2,'lineStyle','--'), hold on
plot(etax1en4,'color',mpg,'lineWidth',2,'lineStyle','-.'), hold on
plot(etax1en5,'color',bla,'lineWidth',2,'lineStyle','-')

hLegend = legend('$\eta$ 0','$\eta$ 0.1','$\eta$ 0.01', '$\eta$ 0.001',...
    '$\eta$ 0.0001','$\eta$ 0.00001');
set(hLegend,'Interpreter','Latex','Location','northeast');
hXLabel = xlabel('$\# steps$');
hYLabel = ylabel('$relative\ error$');
thisFigure
%% s_snr5_etafDiff_etax1e-2
startup;

load('s_snr5_etafDiff_etax1e-2');

figure(2)
plot(etaf1en2,'color',dre,'lineWidth',2,'lineStyle','-'), hold on
plot(etaf1en3,'color',ora,'lineWidth',2,'lineStyle',':'), hold on
plot(etaf1en4,'color',blu,'lineWidth',2,'lineStyle','-.'), hold on
plot(etaf1en5,'color',gra,'lineWidth',2,'lineStyle','--')

hLegend = legend('$\eta$ 0.01', '$\eta$ 0.001',...
    '$\eta$ 0.0001','$\eta$ 0.00001');
set(hLegend,'Interpreter','Latex','Location','northwest');
hXLabel = xlabel('$\# steps$');
hYLabel = ylabel('$relative\ error$');
thisFigure
%% s_snr5_etaf1e-4_etaxDiff(!)
startup;
figPath = '/home/gao/Documents/MPI/thesis/article/figure/fig4thesis';

load('s_snr5_etaf1e-4_etaxDiff');
load('s_snr5_etafDiff_etax1e-2');
load('err_s_snr5_etaf0.01_etax0'); % (!)

frame = 0:100;
figure(1),clf
plot(frame,err.rel,'color',dre,'lineWidth',2,'lineStyle','-'), hold on
plot(frame,etaf1en3,'color',ora,'lineWidth',2,'lineStyle',':'), hold on
plot(frame,etax1en2,'color',blu,'lineWidth',2,'lineStyle','-.'), hold on
plot(frame,etax1en3,'color',gra,'lineWidth',2,'lineStyle','--'), hold on
plot(frame,etax1en5,'color',mpg,'lineWidth',2,'lineStyle','-.')

hLegend = legend('$\eta$ 0','$\eta$ 0.1','$\eta$ 0.01', '$\eta$ 0.001',...
    '$\eta$ 0.00001');
set(hLegend,'Interpreter','Latex','Location','northwest');
hXLabel = xlabel('$\#\text{frames}$', 'Interpreter','Latex');
hYLabel = ylabel('$\text{relative error}$'); 
thisFigure
figName = 's_etaf1e-4_etaxDiff';
figName = strcat(figName,'.tikz');
figname = fullfile(figPath,figName);
printTikz;
%% s_snrInf_etaf0_etax0
startup;

figure(1),clf
plot(0:100,err.rel,'color',dre,'lineWidth',2,'lineStyle','-')

hXLabel = xlabel('$\# frames $');
hYLabel = ylabel('$relative\ error$');
thisFigure
filename = 's_snrInf_etaf0_etax0';
print(gcf, '-deps2c', filename)