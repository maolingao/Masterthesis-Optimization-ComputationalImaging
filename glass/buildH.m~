eps = 1e-30;
% tail = zeros(numel(H.H));
tail = zeros(size(H.H));
for k = 1 : H.i - 1

    den = (H.s(:,k)'*H.y(:,k) );
    tail = tail + (den + eps )\(H.delta(:,k)*(H.s(:,k))' ) + ...
        (den + eps )\(H.s(:,k)*(H.delta(:,k)') ) - ...
        (den^2 + eps )\(H.s(:,k)*(H.delta(:,k)'*H.y(:,k))*(H.s(:,k)') );
end
% H_mtx = tail + kron(H.H,H.H);
H_mtx = tail + H.H;

% den_ref = s_ref'*y_ref;
% H_ref = H_ref + (den_ref)\(delta_ref*s_ref') + ...
%     (den_ref)\(s_ref*delta_ref') - ...
%     ((den_ref)^2)\(s_ref*delta_ref'*y_ref*s_ref');% H_i+1 <-- H_i + (update)

