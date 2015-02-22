function drawComparisonFig(gtImg,cmpImg,counter,optimizerName,cmpTargetName,figPath)
% 
% -------- comparison figure --------
fImg = figure; set(fImg,'visible','off'),
subplot(1,2,1)
imagesc(clip(gtImg,1,0)); 
axis image off; colormap gray;
subplot(1,2,2)
imagesc(clip(cmpImg,1,0)); 
axis image off; colormap gray;
filename = sprintf('Com_%d',counter);
filename = strcat(optimizerName,cmpTargetName,filename);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
close gcf;
%
% -------- comparison figure --------
fImg = figure; set(fImg,'visible','off'),
imagesc(clip(cmpImg,1,0)); 
axis image off; colormap gray;
filename = sprintf('Img_%d',counter);
filename = strcat(optimizerName,cmpTargetName,filename);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
close gcf;
end