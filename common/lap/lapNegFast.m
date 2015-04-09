function Lm_neg = lapNegFast(m)
% calculate negative part of laplacian operator
% L = Lp - Ln
% Ln = [0 1 0;
%       1 0 1;
%       0 1 0];

Lm_neg = zeros(size(m));
shift.s = 1;

shift.p = 'up'; 
m_shift_pad = shiftPad(m,shift); Lm_neg = Lm_neg + m_shift_pad;
shift.p = 'down'; 
m_shift_pad = shiftPad(m,shift); Lm_neg = Lm_neg + m_shift_pad;
shift.p = 'right'; 
m_shift_pad = shiftPad(m,shift); Lm_neg = Lm_neg + m_shift_pad;
shift.p = 'left'; 
m_shift_pad = shiftPad(m,shift); Lm_neg = Lm_neg + m_shift_pad;

end