function Lf = lap(f, domain)
if ~exist('domain', 'var')
  Lf = -4*del2(f);
elseif exist('domain', 'var')
  if domain == '+'
    Lf = 4*f;
  elseif domain == '-'
    Lf = 4*(del2(f) + f);
  else
    print '[laplacian.m]: domain can be either '+' or '-''
  end
end
return 