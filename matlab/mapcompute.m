function Ea = mapcompute(chi, opbasis)
% compute diag of process matrix
[V,D] = eig(chi)
% compute dimension of system
[ind,~] = size(V);
dim = log2(ind);
% unitary matrix that diags chi
U = V';
Ea = zeros(2,2,dim^2);
for i=1:dim^2
    tmp = zeros(dim);
    for j=1:dim^2
        tmp = tmp + U(i,j)*opbasis{j}*sqrt(D(i,i));
    end
    Ea(:,:,i) = tmp;
end
end
