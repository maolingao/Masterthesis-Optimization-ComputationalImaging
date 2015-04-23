% tikztest

figure(1)
x = -2*pi: pi/10 : 2*pi;
y = sin(x);
hData = plot(x,y,'color','r','lineWidth',2,'lineStyle','-.');
hTitle = title('tikz fig');
hXLabel = xlabel('$\theta$');
hYLabel = ylabel('y');
hLegend = legend('sinus');
thisFigure;

figPath = '/home/gao/gitHub/NOCI/fig';
figName = 'tikztest.tikz';
filename = fullfile(figPath,figName);
matlab2tikz(filename,...
	'height','\figheight','width','\figwidth',...
	'parseStrings',false,'dataPath', './fig','relativeDataPath','./fig','imagesAsPng', true,...
	'extraAxisOptions',{'ylabel near ticks'})