function yi = mfbd_load(i, aux)

fnamei = aux.fname(aux.indices(i));
fprintf('processing %s\n', fnamei);
sy = aux.sy; 

% get the next frame
if all(fnamei(end-2:end) == 'fit')
  yi = double(fitsread(fnamei));
else
  yi = double(imread(fnamei));
end
% pick color channel
yi = yi(:,:,aux.colorchannel);

if isfield(aux, 'shiftfile')
  shifts = dlmread(aux.shiftfile);
  shifts = shifts(aux.indices(i),:);
else
  shifts = [i 0 0];
end

if isfield(aux, 'crop')
  % crop a region
  crop = aux.crop;
  if (aux.initx == 4) && (i == 1)
    sf = aux.sf/2;
    yi = yi(crop(1)+(1-sf(1):sy(1)+sf(1)-1)-shifts(2), crop(2)+(1-sf(2):sy(2)+sf(2)-1)-shifts(3), :);
  else
    yi = yi(crop(1)+(1:sy(1))-shifts(2), crop(2)+(1:sy(2))-shifts(3), :);
  end
end


% aux.fname = @(i) sprintf('%s/kopernikus%04d.png', imagepath, i);
% aux.shiftfile = sprintf('%s/alignment.ausw', imagepath);
% aux.indices = [1:2350]; 