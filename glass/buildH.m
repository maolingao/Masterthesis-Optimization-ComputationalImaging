epsll = 1e-9;
% tail = zeros(numel(H.H));
tailM = zeros(size(H.H));

for k = 1 : H.i - 1

    if ~isempty(H.s)
        if size(H.s,2) == size(H.Ginv0,1)
            Ginv = H.Ginv0;
        elseif size(H.s,2) < 2 
            Ginv = 1/((H.s)'*H.s + epsll);
        else
            G = H.s' * H.y;
            %----------------------------%
            % pseudo-inverse
            Ginv = pinv(G);
            %----------------------------%
            % backslash
%             Ginv = G \ eye(size(H.s,2));
        end
        SGinv = H.s * Ginv;
        SGinvDelta = SGinv * (H.delta)';
        tailM = SGinvDelta + SGinvDelta' - SGinv * H.y' * SGinvDelta';
        H.Ginv0 = Ginv;
        %}
    else
        tailM = 0;
    end
end

% *** update form: H = H0 + tail ***
H_mtx = tailM + H.H;

