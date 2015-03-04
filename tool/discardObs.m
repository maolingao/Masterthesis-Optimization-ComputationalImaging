function [S,Y,Delta] =  discardObs(s, y, delta, cutLine)
% discard old observations in s, y, delta
%
%
if cutLine == 0 || isnan(cutLine)
    
    S       =   s;
    Y       =   y;
    Delta   =   delta;
else
    
    S       =   s(:,cutLine:end);
    Y       =   y(:,cutLine:end);
    Delta   =   delta(:,cutLine:end);

end