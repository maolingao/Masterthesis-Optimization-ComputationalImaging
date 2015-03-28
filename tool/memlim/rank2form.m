function [V,U] = rank2form(s,y,delta,option)
% return rank-2 update block component U and V
% which obeys U'*V + U*V' = delta*G^(-1)*s' + s*G^(-1)*delta' - s*G^(-1)*y'*delta*G^(-1)*s'

% compute W*y via handle
if isa(option.Wfun,'function_handle')
    keyboard
    z = option.Wfun(y);
elseif isa(option.Wfun,'char') && strcmp(option.Wfun,'BFGS')
% implicit W=H function:
    z =  s;
elseif isa(option.Wfun,'char') && strcmp(option.Wfun,'GS')
% explicitly W=H0 function:
    if ~isempty(option.data.R) && ~isempty(option.data.D)
        z = option.data.R * option.data.D * (option.data.R' * y);
    else
        z = y;
    end
else
    error('malformed covariance function')
end

G = z' * y;
G = 1/2*(G + G');

Ginv = pinv(G);
V = z * Ginv;
% V = s / G;
U = delta - 1/2 * V * y' * delta;

end