function result = betashape(beta, dim)
result = zeros(dim^4);
%chk = permute(reshape(1:256,[4,4,4,4]), [2,4,1,3]);
beta = permute(beta,[2,4,1,3]);
for j=1:dim^2
    for k=1:dim^2
        for m=1:dim^2
            for n=1:dim^2
                [r,s] = subget(m,j,n,k, dim);
                result(r,s) = beta(m,j,n,k);
            end
        end
    end
end
end

