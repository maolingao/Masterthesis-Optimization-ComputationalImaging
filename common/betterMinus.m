function sumDiff = betterMinus(cellArray1, cellArray2)
% cell array minus
% if cellArrays only have one cell or just matrices, then recover normal matrix substraction
%
if iscell(cellArray1) && iscell(cellArray2)
   % cellarrays
   if numel(unique([length(cellArray1), length(cellArray2)]))
       cellNum = length(cellArray1);
       sumDiff = zeros(size(cellArray1{1}));
   else
       disp('[betterMinus.m] : two cell arrays to minus must have the same length.')
   end
   %
   for i = 1 : cellNum
       sumDiff = sumDiff + ( cellArray1{i} - cellArray2{i} );
   end
elseif ismatrix(cellArray1) && ismatrix(cellArray2)
    % matrices
    sumDiff = cellArray1 - cellArray2;
else
    % 
   disp('[betterMinus.m] : inputs must be the same type, either 2 cellarrays or 2 matrices.')
end
    
end