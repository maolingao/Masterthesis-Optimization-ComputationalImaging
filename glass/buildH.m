epsll = 1e-9;
% tail = zeros(numel(H.H));
tailM = zeros(size(H.H));

for k = 1 : H.i - 1

    if ~isempty(H.s)
        if size(H.s,2) == size(H.Ginv0,1)
%                         keyboard
            Ginv = H.Ginv0;
        elseif size(H.s,2) < 2 
            Ginv = 1/((H.s)'*H.s + epsll);
        else
            keyboard
            Ginv = invGram(obj.Ginv0,obj.s,obj.y);
%             Ginv = ((obj.s)'*obj.y) \ eye(size(obj.s,2));
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

% ########################
% check, update fasion of H
%{
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
H_mtx = tailM + H.H;
% ########################

