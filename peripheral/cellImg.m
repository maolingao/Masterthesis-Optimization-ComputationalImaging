function imgCell = cellImg(mtx,imgSize)
% convert columns in matrix 'mtx' into images of size 'imgSize'
if ~exist('imgSize','var'), imgSize = sqrt(size(mtx,1)).*[1,1]; end

len = size(mtx,2);
imgCell = cell(len,1);
for i = 1 : len
   imgCell{i} = reshape(mtx(:,i),imgSize); 
end

end