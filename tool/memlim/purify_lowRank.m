function [R,D] = purify_lowRank(s,y,delta,MEMLIM,r,d,option)
% wrap up low rank storage limited update

if exist('r','var') && exist('d','var') && ~isempty(r) && ~isempty(d)   % MEMLIM both current observation and R*D*R'
    if ~isempty(s) && ~isempty(y)
        [R,D] = purify_lowRank_2memlims(r,d,s,y,delta,MEMLIM,option);
    else
        R = r; D = d;
    end
elseif ~isempty(s) && ~isempty(y)
    [R,D] = purify_lowRank_obs(s,y,delta,MEMLIM,option);  % MEMLIM current frame H
else
    error('[purify_lowRank.m] : no data input!')

end

end
