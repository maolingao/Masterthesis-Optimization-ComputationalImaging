% ####### figure ######
% save convolution problem setup figure


figPath = option.figPath;
% for debug
figure(14),  set(gcf,'visible','on'), clf
figure(15),  set(gcf,'visible','on'), clf
% for latex
figure(12),  set(gcf,'visible','off'), clf
figure(10),  set(gcf,'visible','off'), clf
figure(11),  set(gcf,'visible','off'), clf
figure(13),  set(gcf,'visible','off'), clf

%----- nature -----
fnature = figure; set(fnature,'visible','off')
imagesc(nature), colormap gray, axis image off, box off
% title('ground truth')
filename = 'natureIm';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
close gcf

%----- kernel -----
fkernel = figure; set(fkernel,'visible','off')
imagesc(kernel), colormap gray, axis image off, box off
% title('kernel')
filename = 'kernelIm';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
close gcf

%----- blur image -----
fblur = figure; set(fblur,'visible','off')
imagesc(convIm), colormap(gray), axis image off, box off
if exist('SNR','var')
    title(sprintf('SNR = %ddB',SNR))
else
%     title(sprintf('blurry image'))
end
filename = 'noiseIm';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
close gcf