function m = mask(x,tau)

if ~exist('tau','var'), tau = 0.3; end

[r_map, numerator, denominator] = rmap(x);

m = heaviside(r_map - tau);

% figure, imagesc(m), colormap gray, axis image off
end