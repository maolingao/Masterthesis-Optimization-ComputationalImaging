function tightSubplot(imgCell, gap, figureName, figPath, counter)
% wrap up tightPlots.m
% subplot images in a cell array beside each other with controllable gap
% between subfigures.
if ~exist('gap','var'), gap = [0 0]; end
if ~exist('figureName','var'), figureName = ''; end
if ~exist('figPath','var'), figPath = '/is/ei/mgao/figure2drag'; end
if ~exist('counter','var'), counter = ''; end

f           =   length(imgCell);
[nrow,ncol] =   factorize(f);
figure;
ha          =   tightPlots(nrow, ncol, 30, [1 1], gap, [0.6 0.7], [0.8 0.3], 'centimeters');
%
maxPixel = 0;
minPixel = 0;
for i = 1 : f                        % max and min pixel value of all imgs
    maxPixel = max([maxPixel; vec(imgCell{i})]);
    minPixel = min([minPixel; vec(imgCell{i})]);
end
for i = 1 : f
    set(gcf,'visible','off');
    axes(ha(i)); imagesc(imgCell{i}), % axis off % image
    caxis manual
    caxis([minPixel maxPixel]);
%     caxis([-0.25 0.25]);
    colormap default
end
hb = colorbar('location','eastoutside');
% keyboard
set(gcf,'visible','off'),
set(ha(1:f), 'fontname', 'Times', 'fontsize', 10)
set(ha(1:f), 'xticklabel', '','xtick',[]);
set(ha(1:f), 'yticklabel', '','ytick',[]);
% axes(ha(1));  title('Title 1');
% axes(ha(2));  title('Title 2');

% print(gcf, 'Example1.png', '-dpng', '-r200', '-opengl');
% print(gcf, 'Example1.eps', '-depsc2', '-painters', '-loose');
% print(gcf, 'Example1.pdf', '-dpdf', '-painters', '-loose');
filename = sprintf('Img_%d',counter);
filename = strcat(figureName,filename);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename, '-painters', '-loose');
close gcf
end
