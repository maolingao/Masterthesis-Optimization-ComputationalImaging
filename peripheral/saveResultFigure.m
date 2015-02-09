% ####### figure ######
% save convolution problem solution and diagnostic figure

figPath = option.figPath;

%% plot
% for debug
fclk = figure(14); set(fclk,'visible','on'),
subplot(121)
hTitle = title(shape);
% legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
hYLabel = ylabel ('$\|Fx - y\| / pixel$','Interpreter','Latex');
hXLabel = xlabel('$time/sec$','Interpreter','Latex');
thisFigure;   
subplot(122)
hTitle = title(shape);
% legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
hYLabel = ylabel ('$\|x - \hat{x}\| / \|x\|$','Interpreter','Latex') ;
hXLabel = xlabel('$time/sec$','Interpreter','Latex');
thisFigure;   
%
fstp = figure(15); set(fstp,'visible','on'),
subplot(121)
hTitle = title(shape);
% legend('my pncg','my cg','my lucy','my gaussian');legend boxoff
hYLabel = ylabel ('$\|Fx - y\| / pixel$','Interpreter','Latex');
hXLabel = xlabel('$\#steps$','Interpreter','Latex');
thisFigure;   
subplot(122)
hTitle = title(shape);
% legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
hYLabel = ylabel ('$\|x - \hat{x}\| / \|x\|$','Interpreter','Latex');
hXLabel = xlabel('$\#steps$','Interpreter','Latex');
thisFigure;   
%%
if isfield(option,'mode')
    switch option.mode
        case 'compPerFrame'
            counter = j;
        case 'compAllFrame'
            counter = inf;
        otherwise
            display('in [mbd.m]: option.mode can be either "compPerFrame" or "compAllFrame"');
    end
else
    counter = inf;
end
    

%% plot
% for latex
f10 = figure(10); set(f10,'visible','off'),
hTitle = title(shape);
% legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
hYLabel = ylabel ('$\|Fx - y\| / pixel$','Interpreter','Latex');
hXLabel = xlabel('$time/sec$','Interpreter','Latex');
thisFigure;   
filename = sprintf('compClk_residual_frame_%d',counter);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
close gcf;
%
f12 = figure(12); set(f12,'visible','off'),
hTitle = title(shape);
% legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
hYLabel = ylabel ('$\|x - \hat{x}\| / \|x\|$', 'Interpreter','Latex'); % set(gca,'yaxislocation','right');
% p=get(ylh,'position'); set(gca,'yaxislocation','left');  set(ylh,'position',p);
hXLabel = xlabel('$time/sec$','Interpreter','Latex');
thisFigure;   
filename = sprintf('compClk_relErr_frame_%d',counter);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
close gcf;
%  
f11 = figure(11); set(f11,'visible','off'),
hTitle = title(shape);
% legend('my pncg','my cg','my lucy','my gaussian');legend boxoff
hYLabel = ylabel ('$\|Fx - y\| / pixel$','Interpreter','Latex');
hXLabel = xlabel('$\#steps$','Interpreter','Latex');
thisFigure;   
filename = sprintf('compStp_residual_frame_%d',counter);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
close gcf;
%
f13 = figure(13); set(f13,'visible','off'),
hTitle = title(shape);
% legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
hYLabel =  ylabel ('$\|x - \hat{x}\| / \|x\|$','Interpreter','Latex');  %set(gca,'yaxislocation','right');
% p=get(ylh,'position'); set(gca,'yaxislocation','left');  set(ylh,'position',p);
hXLabel = xlabel('$\#steps$','Interpreter','Latex');
thisFigure;   
filename = sprintf('compStp_relErr_frame_%d',counter);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
close gcf;
% close all hidden figures
% close all hidden