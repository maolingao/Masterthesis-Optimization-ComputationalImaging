function [L,R] = reduce(LH,RH, MEMLIM)
% Memory limited update, keep only dominant information, obtained from
% Maren
    zlr = LH*RH';
    if norm(zlr) == 0 % if low rank term is zero
        L = LH;
        R = RH;
    else
        [QL, SL, VL] = svd(LH, 'econ');
        [QR, SR, VR] = svd(RH, 'econ');
        [q, D, v]    = svd(SL*VL'*VR*SR', 'econ');

        Q   = QL*q;                    % compute svd of low rank term: 
        V   = QR*v;                    % [Q, S, V] = svd(L*R', 'econ')
        d   = diag(D);

        idx = find(d > 1e-6);          % singular values greater than threshold
        keyboard
        if isempty(idx)
            L = LH;
            R = RH;
            return;
        end
        if length(idx) > MEMLIM
            idx = idx(1:MEMLIM);       % cut off at memory limit
        end
        Q   = Q(:, idx);               % reduce rank
        V   = V(:, idx); 
        d   = d(idx);
        L = Q*diag(d);                 % make instance of reduced matrix 
        R = V;        
    end
end
