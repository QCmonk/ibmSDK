function nums = z2u(op)
cvx_begin quiet
    cvx_precision best
    variable nums(1,3)
    minimise( norm(u3(nums) - op))
    subject to
        0 <= nums <= 1;
cvx_end
end