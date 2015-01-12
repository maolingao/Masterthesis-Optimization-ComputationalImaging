% ####### figure ######
% save convolution problem setup figure


figPath = option.figPath;
figure(12), clf
figure(10), clf
figure(11), clf
figure(13), clf

%----- nature -----
fnature = figure; set(fnature,'visible','off')
imagesc(nature), colormap gray, axis image off, box off
title('ground truth')
filename = 'natureIm';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

%----- kernel -----
fkernel = figure; set(fkernel,'visible','off')
imagesc(kernel), colormap gray, axis image off, box off
title('kernel')
filename = 'kernelIm';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

%----- blur image -----
fblur = figure; set(fblur,'visible','off')
imagesc(convIm), colormap(gray), axis image off, box off
if exist('SNR','var')
    title(sprintf('SNR = %ddB',SNR))
else
    title(sprintf('blurry image'))
end
filename = 'noiseIm';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)