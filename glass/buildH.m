eps = 1e-30;
% tail = zeros(numel(H.H));
tail = zeros(size(H.H));
for k = 1 : H.i - 1

    den = (H.s(:,k)'*H.y(:,k) );
%     display(sprintf('den: %d',den))
    tail = tail + (den + eps )\(H.delta(:,k)*(H.s(:,k))' ) + ...
        (den + eps )\(H.s(:,k)*(H.delta(:,k)') ) - ...
        (den^2 + eps )\(H.s(:,k)*(H.delta(:,k)'*H.y(:,k))*(H.s(:,k)') );
end

% ########################
% check, update fasion of H
%{%
% *** update form: H = H0 + l*r' + tail ***
if ~isempty(H.l)
%     keyboard
    lrt = H.l * H.r';
    H_mtx = H.H + lrt + tail;        % output = I + l*r' + tail
else
    % H_mtx = tail + kron(H.H,H.H);
    H_mtx =  H.H + tail;            % output = I + tail
end
%}
% *** update form: H = H0 + tail ***
% H_mtx = tail + H.H;
% ########################

