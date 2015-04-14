function Lf = lap(f, domain)
% wrap-up lapFast.m

if ~exist('domain', 'var')
%   Lf = -4*del2(f);
  Lf = lapFast(f);
elseif exist('domain', 'var')
  if domain == '+'
    Lf = 4*f;
  elseif domain == '-'
%     Lf = 4*(del2(f) + f);
    Lf = lapNegFast(f);
  else
    print '[lap.m]: domain can be either '+' or '-''
  end
end
return 