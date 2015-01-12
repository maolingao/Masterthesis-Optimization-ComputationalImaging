% ####### figure ######
% save convolution problem solution and diagnostic figure

figPath = option.figPath;

%% plot
% for debug
fclk = figure(14); set(fclk,'visible','on'),
subplot(121)
title(shape);
legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
ylabel('$\|Fx - y\| / pixel$')
xlabel('$time/sec$')
subplot(122)
title(shape);
legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
ylabel('$\|x - \hat{x}\| / \|x\|$'); 
xlabel('$time/sec$')
%
fstp = figure(15); set(fstp,'visible','on'),
subplot(121)
title(shape);
legend('my pncg','my cg','my lucy','my gaussian');legend boxoff
ylabel('$\|Fx - y\| / pixel$')
xlabel('$\#steps$')
subplot(122)
title(shape);
legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
ylabel('$\|x - \hat{x}\| / \|x\|$');
xlabel('$\#steps$')
% for latex
f10 = figure(10); set(f10,'visible','off'),
title(shape);
legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
ylabel('$\|Fx - y\| / pixel$')
xlabel('$time/sec$')
filename = ['compClk_residual'];
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
%
f12 = figure(12); set(f12,'visible','off'),
title(shape);
legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
ylh = ylabel('$\|x - \hat{x}\| / \|x\|$'); % set(gca,'yaxislocation','right');
% p=get(ylh,'position'); set(gca,'yaxislocation','left');  set(ylh,'position',p);
xlabel('$time/sec$')
filename = ['compClk_relErr'];
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
%  
f11 = figure(11); set(f11,'visible','off'),
title(shape);
legend('my pncg','my cg','my lucy','my gaussian');legend boxoff
ylabel('$\|Fx - y\| / pixel$')
xlabel('$\#steps$')
filename = ['compStp_residual'];
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
%
f13 = figure(13); set(f13,'visible','off'),
title(shape);
legend('my pncg','my cg','my lucy','my gaussian'); legend boxoff
ylh = ylabel('$\|x - \hat{x}\| / \|x\|$');  %set(gca,'yaxislocation','right');
% p=get(ylh,'position'); set(gca,'yaxislocation','left');  set(ylh,'position',p);
xlabel('$\#steps$')
filename = ['compStp_relErr'];
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
% close all hidden figures
% close all hidden