function [R,D] = purify_lowRank_obs(s,y,delta,MEMLIM,option)
% low rank matrix storage limited purify, based on observations within
% linear problems


% compute W*y via handle
if isa(option.Wfun,'function_handle')
    z = option.Wfun(obj.y);
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

Z = 1/sqrt(2) * [V - U, V + U];
i = ones(size(V,2),1);
e = [-i; i];
% E = diag(-i,i); Ztilde = E * Z';

[R,D] = Eig_LowRankSymmetricRealMatrix(Z,e,MEMLIM);



end