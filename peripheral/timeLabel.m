function time_accu = timeLabel(time_accu, time_single)
% assemble overall time label from single time  label
if time_accu == 0
    time_accu = [];
end
if isempty(time_accu)
    time_accu = time_single;
else
    time_single = time_accu(end) + time_single;
    time_accu = [time_accu, time_single];
end

end