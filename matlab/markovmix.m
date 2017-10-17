function lambda = markovmix(basestate, markov)
% yep
beep off
% get size of spanning set
basenum = length(basestate);
% get dimension of states
[basedim,~] = size(markov); 
% get number of timesteps
steps = log2(basedim)/2;
%reshape basestates as required
spanset = zeros(basenum, basedim^2);
for i=1:basenum
    spanset(i,:) = reshape(basestate{i}, 1, basedim^2);
end

% begin optimisation 
cvx_begin sdp
    cvx_precision best
    variable lambda(1, basenum)
    minimise(norm( reshape(lambda*spanset, basedim, basedim) - markov ))
    subject to
        % Choi state is positive
        reshape(lambda*spanset, basedim, basedim) >= 0;
        % Choi state is a valid density matrix
        trace(reshape(lambda*spanset, basedim, basedim)) == 1;
cvx_end

end