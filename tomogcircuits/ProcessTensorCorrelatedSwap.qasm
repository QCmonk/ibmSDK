OPENQASM 2.0;

include "qelib1.inc";
qreg q[5];
creg c[5];



h q[0];
h q[3];
cx q[0], q[1];
cx q[3], q[2];

s q[0];

cx q[0],q[2];
h q[0];
h q[2];
cx q[0], q[2];
h q[0];
h q[2];
cx q[0], q[2];

cx q[2],q[0];

t q[0];
