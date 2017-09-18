function mu = dnorm2(N,A)
beep off
% Computes the diamond norm betwen the two quantum channels N,A in B form. 
[dim, ~] = size(N);
cvx_begin sdp
    cvx_precision best
    variable mu 
    variable ZAB(dim,dim) complex
    minimize(mu)
    subject to:
        ptrace(ZAB, 2, [dim/2,dim/2]) <= mu*eye(dim/2);
        ZAB >= N - A;
        ZAB >= 0;
cvx_end
end

