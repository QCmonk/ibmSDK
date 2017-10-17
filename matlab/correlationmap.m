function corrset = correlationmap(ptensor, ptensormarkov, qre)
beep off

% dimension of system
[dim,~] = size(ptensor);

% compute correlation map
cvx_begin sdp
    % fucking fine CVX, I'll do it your way
    variable G(dim^2,dim^2) hermitian semidefinite
    maximise(sqrt(-(quantum_rel_entr(ptensormarkov, reshape(choi_liou_involution(G)*reshape(ptensor, dim^2,1),dim,dim)) - qre)))
    subject to 
        % new state is positive definite
        0.5*(reshape(choi_liou_involution(G)*reshape(ptensor, dim^2,1),dim,dim) + reshape(choi_liou_involution(G)*reshape(ptensor, dim^2,1),dim,dim)') >= 0;
        % new state has trace 1
        trace(0.5*(reshape(choi_liou_involution(G)*reshape(ptensor, dim^2,1),dim,dim) + reshape(choi_liou_involution(G)*reshape(ptensor, dim^2,1),dim,dim)')) == 1;
        % choi state is trace preserving
        ptrace(G, [1], [dim,dim]) == eye(dim);
cvx_end
corrset = G;
end
