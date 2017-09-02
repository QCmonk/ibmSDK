function chi = processtomography(reconstruct, densitybasis, operatorbasis)
beep off
% compute dimension of system
[dim, ~] = size(reconstruct{1});
%--------------------------------------------------------------------------
%-------------------------- STAGE 1 OPTIMISATION---------------------------
%--------------------------------------------------------------------------
% Compute Lambda matrix

A = zeros(dim^2);
b = zeros(dim^2);
for i=1:dim^2
    A(i,:) = reshape(densitybasis{i}.', [1,dim^2]);
    b(i,:) = reshape(reconstruct{i}, [1,dim^2]);
end

% Lam is a two dimensional array with columns indexing input basis states 
% and rows indexing terms in the sum: sum_{k} lam_{jk}*rho_k
cvx_begin quiet
        cvx_precision best
        variable lam(dim^2, dim^2) complex
        minimise( norm(lam.'*A - b))
cvx_end


%--------------------------------------------------------------------------
%-------------------------- STAGE 2 OPTIMISATION---------------------------
%--------------------------------------------------------------------------
% compute beta matrix
    
% flatten operators
op = zeros(dim^2);
rho = zeros(dim^2);
for i=1:dim^2
   op(:,i) = reshape(operatorbasis{i}.', [1,dim^2]);
   rho(i,:) = reshape(densitybasis{i}.', [1,dim^2]);
end

% construct constrained basis set acting on basis states baserho(m,j,n, rho)
baserhobase = zeros(dim^4);

for m=1:dim^2
    for j=1:dim^2
        for n=1:dim^2
            [r,s] = subget(m,j,n,1:dim^2,dim);
            baserhobase(s,r) = reshape((operatorbasis{m}*densitybasis{j}*operatorbasis{n}').', [1,dim^2]);
        end
    end
end


% perform state estimation for beta matrix - are there constraints I can
% use?
betab = zeros(dim^4);
for j=1:dim^2
    ind = (j-1)*dim^2 +[1:dim^2];
    cvx_begin quiet
        cvx_precision high
        variable bbeta(dim^4, dim^2) complex
        minimise( norm(bbeta*rho - baserhobase(:,ind)))
    cvx_end
    betab(ind,:) = bbeta.';
end

%--------------------------------------------------------------------------
%-------------------------- STAGE 3 OPTIMISATION---------------------------
%--------------------------------------------------------------------------
% compute the Chi matrix

% compute base products for use in constraints
basecont = [];
for m=1:dim^2
    for n=1:dim^2
        % pre allocate this in future version, will need to sort out indexing
        basecont = vertcat(basecont, operatorbasis{n}'*operatorbasis{m});
    end
end

cvx_begin
    cvx_precision best
    variable chi(dim^2,dim^2) hermitian semidefinite
    minimise( norm( (pinv(full(betab))*reshape(lam, [dim^4,1])) - reshape(chi, [dim^4,1])))
    subject to
        % trace preserving (or not) constraint (pain in the neck to do in one step...CVXXXXX! *shakes fist*)
        real(repmat(eye(dim),1,dim^4)*(interleave(repmat(reshape(chi, [dim^4,1]),[1,dim]),dim).*basecont)) <= eye(dim);
        imag(repmat(eye(dim),1,dim^4)*(interleave(repmat(reshape(chi, [dim^4,1]),[1,dim]),dim).*basecont)) == 0;
        % trace relation between chi and P
        trace(chi) == trace(repmat(eye(dim),1,dim^4)*(interleave(repmat(reshape(chi, [dim^4,1]),[1,dim]),dim).*basecont))/dim;
cvx_end

% python is not a fan of sparse matrices
chi = full(chi);
end