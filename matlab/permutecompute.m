function permutemat = permutecompute(ptensororig, ptensortarg)
% compute dimension of system
[dim, ~] = size(ptensororig);

cvx_begin
    cvx_precision best
    variable permutemat(dim,dim)
    minimise(norm( permutemat*ptensororig - ptensortarg ))
cvx_end
end