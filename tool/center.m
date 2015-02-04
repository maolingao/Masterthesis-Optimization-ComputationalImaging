function g = center(f, support);
% usage: shift image to center
%
% -- input --
% f         : image to be shifted
% support   : type of center of mass, 1-geometrical, 0-weighted geometrical
% -- output --
% g         : shifted image
%
%
if ~exist('support','var') || isempty('support'), support = 0; end
cm = centerofmass(f, support);
sf = size(f);
g = shift(f, round(sf/2 - cm));
return
