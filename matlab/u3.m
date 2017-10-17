function U = u3(nums)
theta = nums(1);
phi = nums(2);
lambda = nums(3);

U = [cos(theta/2), -exp(1j*lambda)*sin(theta/2);
     exp(1j*phi)*sin(theta/2), exp(1j*lambda + 1j*phi)*cos(theta/2)];

end

