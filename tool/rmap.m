function [r_map, numerator, denominator] = rmap(x)

[gu,gv] = gradient(x);
hs = 30;
h = ones(hs);
sx = size(x);
Nh = conv2MatOp(h, sx, 'same');
lambda = 0.5;

numerator = sqrt((Nh*gu).^2 + (Nh*gv).^2 );
denominator = Nh*sqrt(gu.^2 + gv.^2);
denominator= denominator + lambda;

r_map = numerator./(denominator+eps);
r_map = r_map./max(vec(r_map));

% figure
% subplot(1,3,1), imagesc(numerator), colormap gray, title('num'), axis image off
% subplot(1,3,2), imagesc(denominator), colormap gray, title('den'), axis image off
% subplot(1,3,3), imagesc(r_map), colormap gray, title('r-map'), axis image off

end