% investigates the new possible modelling method

% define the giant fucking basis set for single qubit unitaries
z = zeros(2,2,10);
z(:,:,1) = eye(2);
for i=2:4
   z(:,:,i) = (eye(2) + 1j*pauli(i-1))/sqrt(2);
   z(:,:,i+3) = (eye(2) - 1j*pauli(i-1))/sqrt(2);
end
z(:,:,8) = (eye(2) + pauli(1)*1j/sqrt(2) + pauli(2)*1j/sqrt(2))/sqrt(2);
z(:,:,9) = (eye(2) + pauli(1)*1j/sqrt(2) + pauli(3)*1j/sqrt(2))/sqrt(2);
z(:,:,10) = (eye(2) + pauli(2)*1j/sqrt(2) + pauli(3)*1j/sqrt(2))/sqrt(2);

(eye(2) + pauli(1)*1j/sqrt(2) + pauli(2)*1j/sqrt(2))/sqrt(2)
(eye(2) + pauli(1)*1j/sqrt(2) + pauli(3)*1j/sqrt(2))/sqrt(2)
(eye(2) + pauli(2)*1j/sqrt(2) + pauli(3)*1j/sqrt(2))/sqrt(2)

% define some unitary operator
U = unitary3(0.5*pi, 0*pi, 0.5*pi);

% define choi state equivalent of map:
choiU = kraus2choi(U);
choiU = choiU/trace(choiU);

% compute the choi equivalents of the spanning set
zChoi = zeros(4,4,10);
for i=1:10
    zChoi(:,:,i) = kraus2choi(z(:,:,i));
end

% flatten spanning choi
zChoiflat = zeros(10, 16);
for j = 1:10
   zChoiflat(j, :) = reshape(zChoi(:,:,j), 1,16);
end

cvx_begin sdp
    cvx_precision best
    variable lambda(1, 10)
    minimise(norm(  reshape(lambda*zChoiflat ,4,4) - choiU ))
    subject to
        % positive
        reshape(lambda*zChoiflat ,4,4) >= 0;
        % trace 1
        trace(reshape(lambda*zChoiflat ,4,4)) == 1;
        % unital
        ptrace(reshape(lambda*zChoiflat ,4,4), [2], [2,2]) == eye(2)/2;
cvx_end

lambda
