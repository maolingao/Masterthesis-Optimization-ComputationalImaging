function [R,D] = purify_lowRank(s,y,delta,MEMLIM,r,d,option)
% wrap up low rank storage limited update

[R1,D1] = purify_lowRank_obs(s,y,delta,MEMLIM,option);  % MEMLIM current frame 

if exist('r','var') && exist('d','var') && ~isempty(r) && ~isempty(d)                % aggregate 2 MEMLIMs: current frame MEMLIM and all previous MEMLIM  
    [R,D] = purify_lowRank_2memlims(r,d,R1,D1,MEMLIM);
else
    R = R1;
    D = D1;
    return
end

end