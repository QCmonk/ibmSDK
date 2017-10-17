function U = unitary3(theta, phi, lambda)
U = [cos(theta/2), -exp(1j*lambda)*sin(theta/2);
     exp(1j*phi)*sin(theta/2), exp(1j*lambda + 1j*phi)*cos(theta/2)];
end

