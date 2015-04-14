function m_shift_pad = shiftPad(m,shift)
% shift matrix m and pad empty area with zeros
% shift.p : shift direction : up, down, right, left
% shift.s : shift size : how many rows or columns to be shifted. should be
% smaller than the size(m)


if ~exist('shift','var'); 
    m_shift_pad = m; disp('no shift and pad.'); return
end

ms = size(m); ms_row = ms(1); ms_col = ms(2);
m_shift_pad = zeros(ms);

switch shift.p % shift direction
    case 'up'
        m_shift_pad(1 : ms_row-shift.s, :) = m(1+shift.s : end, :);
    case 'down'
        m_shift_pad(1+shift.s : end, :) = m(1 : ms_row-shift.s, :);
    case 'right'
        m_shift_pad(:, 1+shift.s : end) = m(:, 1 : ms_col-shift.s);
    case 'left'
        m_shift_pad(:, 1 : ms_col-shift.s) = m(:, 1+shift.s : end);
    otherwise
        error('malformed shfit direction.')
end

end