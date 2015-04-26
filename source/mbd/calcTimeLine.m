function timeLine = calcTimeLine(dataK)

setnum = length(dataK);

timeLine = zeros(setnum,1);
timeLine(1) = dataK{1}.time(1);

for k = 1 : setnum
    timeLine(k+1) = timeLine(k) + dataK{k}.time(end);

end


end