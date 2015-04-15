function yi = mfbd_load_wrap(aux)
if nargin < 1
    
    aux.shift = [0,0];
    aux.start = [1346, 1:9];
end

imagepath = '/is/ei/mgao/Documents/thesis/Astro/real_data/sequence';
figPath = '/is/ei/mgao/Documents/thesis/Astro/real_data/copernicus';

aux.fname = @(i) sprintf('%s/kopernikus%04d.png', imagepath, i);
aux.shiftfile = sprintf('%s/alignment.ausw', imagepath);
aux.indices = [1:2350]; 
% aux.shift = [0,0];
aux.crop = [390, 485] - aux.shift;
% aux.crop = [470, 770] - aux.shift;
aux.initx = 0;
aux.sy = [70,70] + 2 * aux.shift; 
aux.colorchannel = 1;

idx = aux.start(2:end);
yi = cell(length(aux.start),1);
yi{1} = mfbd_load(aux.start(1), aux);
for i = idx % aux.indices
yi{i+1} = mfbd_load(i, aux);

% -------- single figure --------
% fImg = figure; set(fImg,'visible','off'),
% imagesc(clip(yi{i},inf,0)); 
% axis image off; colormap gray;
% filename = sprintf('copernicus_%d',i);
% filename = fullfile(figPath,filename);
% print(gcf, '-dpng','-r0', filename)
% close gcf;
end
% figure, imagesc(yi),colormap gray

end