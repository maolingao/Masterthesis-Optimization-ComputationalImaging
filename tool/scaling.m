function [image_scaled] = scaling(image, factor)
% scale image total energy level

if ~exist('factor','var')
    factor = 1e-3;
end

image_scaled = cellfun(@(x) x*factor, image, 'UniformOutput', false);

end