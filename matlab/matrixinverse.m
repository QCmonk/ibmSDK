function inverse = matrixinverse(matrix)
% compute dimension of system
[dim, ~] = size(matrix);

cvx_begin quiet
    cvx_precision best
    variable inverse(dim, dim) complex
    minimise( norm(matrix*inverse - eye(dim)))
cvx_end

end